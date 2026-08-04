# =============================================================================
# compareReplicatesAndMultiIntegrations
# =============================================================================
# Visualization/QC script, not a numbered pipeline stage -- run by hand against stage
# 07/10's output whenever you want this report, same convention as
# plotSplicingEfficiencies_uniform.R and computeIntegrationStats.R.
#
# Verbatim port of igoru's ~/R/trip/compareReplicatesAndMultiIntegrations_uniform.R, no
# logic changes beyond what's already in that source file. QC report per construct: how
# much do independent barcodes landing in the same gene / same intron / same exact position
# agree with each other, compared to random pairs? (expression distance boxplots with
# Wilcoxon comparisons), how many barcodes share a gene/intron/position (histograms), how
# their exp.log10 *levels* (not pairwise distances) compare to a random sample across the
# same four groupings, how much of exp.log10's variance each grouping explains (one-way
# ANOVA eta-squared), and which libraries express which barcodes (UpSet plot, capped to the
# top 40 intersections). Same reads-annotated_tables.txt-plus-combined_counts.txt-per-
# expName merge pattern as comparePairwiseSharedPositions_uniform.R/
# annotateTRIP_positions_genuine_barcodes.R/annotateTRIPGeneBodies_genuine_barcodes.R, but
# the actual comparisons (gene/intron/position multi-integration agreement) are unique to
# this script.
#
# Refreshed 2026-08-04 to pick up three fixes made upstream in igoru's ~/R/trip/ copy
# (previous port here predated all three):
#   - "Shared position" is now strand-specific (Chr:Pos:Strand, not just Chr:Pos) -- two
#     integrations at the same coordinate on opposite strands are distinct events, not a
#     shared position.
#   - Memory-blowout fix #1: the "Random pairs" reference distribution used to sample 10,000
#     barcodes and then take dist() -- every pairwise distance among them, i.e. C(10000,2)
#     ~= 50 million pairs, not the ~10,000 pairs the variable name/label implied. Now
#     samples 10,000 pairs directly.
#   - Memory-blowout fix #2: the UpSet plot's min_size=0 doesn't cap anything -- it renders
#     every *observed* set-combination (up to 2^n_libraries). With enough libraries and
#     barcodes, most of the possible combinations are observed at least once, which is what
#     was actually exhausting memory (confirmed via single-map, single-fix-at-a-time
#     bisection against real SON_combined data -- present in the unmodified pre-strand-fix
#     script too, not caused by the strand-key change). Now capped via n_intersections=40
#     (top 40 largest intersections, standard UpSet convention).
#
# Usage
#   Rscript compareReplicatesAndMultiIntegrations_uniform.R <path> <mapName1,mapName2,...> <expListFile> <outFile.pdf> [titleExtra]
#
# Real example (as actually run for IntronTRIP_combined):
#   bsub -R rusage[mem=30000] Rscript ~/R/trip/compareReplicatesAndMultiIntegrations_uniform.R \
#     ./ A1,B1,C1,D1,E1,B3-6,B3-9,H1 expression.exp_list.txt \
#     plots/replicates_integrations_IntronTRIP_expression.pdf
#
# expListFile: mapName<TAB>expName, one row per (construct, expression-readout) pair to
# include -- like comparePairwiseSharedPositions_uniform.R, real production expListFiles are
# often hand-curated and NOT 1:1 with ### COMBOS (e.g. only the "Total" readout per rep, as
# in the real expression.exp_list.txt this script was designed around). See
# examples/exp_lists/ for real ones, or build one yourself from ### COMBOS' mapName+expSample
# columns.
#
# Inputs (per mapName)   <mapName>.txt.annotated_tables.txt              (stage 10)
# Inputs (per expName)   <expName>.combined_counts.txt                   (stage 07)
# Output                 one multi-page PDF (outFile), 5 pages per mapName: barcodes-per-
#                        gene/intron/position histograms, same-gene/intron/position vs
#                        random-pairs expression-distance boxplots, exp.log10-level
#                        distributions across the same groupings, %-variance-in-exp.log10-
#                        explained by each grouping, and an UpSet plot of which libraries
#                        express which barcodes.
# =============================================================================

# NOTE: on the wexac cluster, `module load R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1`
# crashes the UpSet plot at the end of this script (ggplot2 4's stricter theme validation
# rejects something ComplexUpset/patchwork set internally: "the axis.title.x theme element
# must be a <element_text> object"), confirmed unrelated to this script's own additions
# (fails identically with the unmodified original). `module load R/4.3.1-foss-2022b`
# (older ggplot2) does not have this problem -- use that module to run this script.
library(LSD)
library(ggplot2)
library(ggpubr)
library(ComplexUpset)

# Two barcodes landing at the exact same genomic position are more often PCR/sequencing
# error variants of one true barcode than genuinely independent integrations. Greedily
# keep only barcodes that are mutually >=minDist apart (Hamming, equal-length strings
# only -- a length mismatch is treated as "distant", not merged) within a position group.
dedupBarcodesHamming = function(df, minDist = 3) {
  n = nrow(df)
  if (n <= 1) return(df)
  bcs = strsplit(df$Barcode, "")
  keep = rep(TRUE, n)
  for (i in 1:(n - 1)) {
    if (!keep[i]) next
    for (j in (i + 1):n) {
      if (!keep[j]) next
      if (length(bcs[[i]]) != length(bcs[[j]])) next
      if (sum(bcs[[i]] != bcs[[j]]) < minDist) keep[j] = FALSE
    }
  }
  df[keep, ]
}

# Mean pairwise |av_i - av_j| distance within a gene, restricted to pairs that are NOT
# in the same intron (so "Same gene" and "Same intron" are disjoint comparisons). Two
# rows are treated as "different introns" whenever at least one has NA MANE_intron
# (i.e. not intronic at all) or their non-NA intron IDs differ -- only an exact match of
# two non-NA intron IDs counts as "same intron" and gets excluded. Returns NA if the gene
# has no valid different-intron pair (e.g. all its integrations fall in one intron).
meanDistDiffIntron = function(x) {
  n = nrow(x)
  if (n < 2) return(NA)
  cmb = combn(n, 2)
  dists = c()
  for (k in 1:ncol(cmb)) {
    i = cmb[1, k]; j = cmb[2, k]
    sameIntron = !is.na(x$MANE_intron[i]) && !is.na(x$MANE_intron[j]) && x$MANE_intron[i] == x$MANE_intron[j]
    if (!sameIntron) dists = c(dists, abs(x$av[i] - x$av[j]))
  }
  if (length(dists) == 0) return(NA)
  mean(dists)
}

# One-way ANOVA of a group's av (exp.log10) values against group identity (gene/intron/
# position) -- eta-squared (SS_between/SS_total) is the standard "% of variance explained"
# figure for this, and the same ANOVA's F-test p-value is the appropriate significance test
# for whether that grouping explains a real (non-random) share of the variance. groupList is
# one of spGeneMulti/spIntronMulti/spCompositePos (a list of data.frames, each one gene's/
# intron's/position's barcodes, each with an `av` column) -- every barcode across every
# group is pooled into one long-format table first.
varianceExplained = function(groupList) {
  if (length(groupList) < 2) return(c(pct = NA, p = NA, n = 0, nGroups = length(groupList)))
  df = do.call(rbind, lapply(names(groupList), function(k) {
    data.frame(group = k, value = groupList[[k]]$av)
  }))
  df = df[!is.na(df$value), ]
  fit = tryCatch(aov(value ~ factor(group), data = df), error = function(e) NULL)
  if (is.null(fit)) return(c(pct = NA, p = NA, n = nrow(df), nGroups = length(unique(df$group))))
  tab = summary(fit)[[1]]
  ssBetween = tab["factor(group)", "Sum Sq"]
  ssTotal = sum(tab[, "Sum Sq"])
  c(pct = 100 * ssBetween / ssTotal, p = tab["factor(group)", "Pr(>F)"],
    n = nrow(df), nGroups = length(unique(df$group)))
}

path = "/home/projects/ulitsky/shared/data/human/TRIP/SON_combined///"

mapNames = c("SON_A1","SON_E1","SON_B3_9")
expListFile = "spliced_unspliced.exp_list.txt"
outFile = "plots/replicates_integrations_SON_splicing.pdf"
titleExtra = ""
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  mapNames = strsplit(myArgs[2], ",")[[1]]
  expListFile = myArgs[3]
  outFile = myArgs[4]
  if (length(myArgs)>=5) titleExtra = myArgs[5]
}

expList = read.table(paste0(path,expListFile),sep="\t")

pdf(paste0(path,outFile), width = 12)

for (ind in 1:length(mapNames))
{
  mapName=mapNames[ind]
  print(mapName)
  annotTab=read.table(paste0(path,mapName,".txt.annotated_tables.txt"),sep="\t",comment="")
  dim(annotTab)
  tabs = list()
  allMerged = annotTab
  expNames = expList[expList$V1==mapName,"V2"]
  # Per-library flag: was exp.log10 already a column in the file (e.g. a bounded
  # splicing-ratio metric we don't control the semantics of), or did we compute it
  # ourselves as log10(1+exp)? Only in the latter case is it safe to later fill NAs
  # with a log10(1+exp)-style fallback -- see the post-merge loop below.
  explog10Preexisting = setNames(rep(FALSE, length(expNames)), expNames)
  for (nameExp in expNames)
  {
    countsTab =  read.table(paste0(path, nameExp, ".combined_counts.txt"), sep = "\t", quote = "")
    countsTab$Barcode = rownames(countsTab)
    annotTab$Barcode = rownames(annotTab)
    # If the file already has an exp.log10 column, trust it as-is rather than overwriting
    # it -- only compute it when it's missing.
    explog10Preexisting[nameExp] = "exp.log10" %in% names(countsTab)
    if (!explog10Preexisting[nameExp]) {
      countsTab$exp.log10 = log(1+countsTab$exp,10)
    }
    countsTab$exp.norm = log((1 + countsTab$exp) / (1 + countsTab$norm), 2)
    # Rename by column name, not position: a fixed index like c(1,2,6,8:9) silently
    # renamed the wrong column (or left Barcode/Chr/Pos/Strand unprotected) whenever the
    # file's column count/order differs from the "usual" 6 -- e.g. exactly what happens
    # when exp.log10 already exists as an extra column above.
    for (col in c("norm","exp","map","exp.log10","exp.norm","exp_spliced","exp_total","exp_unspliced")) {
      if (col %in% names(countsTab)) names(countsTab)[names(countsTab)==col] = paste0(nameExp,".",col)
    }
    allMerged = merge(allMerged,countsTab,by=c("Barcode","Chr","Pos","Strand"),all.x=TRUE)
    allMerged$Row.names=NULL
    tabs[[nameExp]]=countsTab
  }
  allMerged$av=apply(allMerged[,grep("exp.log10$",names(allMerged))],1,mean,na.rm=T)

  usable = subset(allMerged,!is.na(allMerged$av))

  # Strand-specific: two integrations at the same Chr:Pos but opposite strands are distinct
  # events, not a shared position, so Strand is part of the grouping key.
  usable$CompositePos = paste0(usable$Chr,":",usable$Pos,":",usable$Strand)

  print(paste(length(unique(usable$MANE_all))," unique genes"))
  print(paste(length(unique(usable$MANE_intron))," unique introns"))

  spGene = split(usable,usable$MANE_all)
  spGeneMulti=spGene[sapply(spGene,function(x){length(unique(x$Pos))>1 & length(unique(x$Pos))==nrow(x)})]
  print(paste(length(spGeneMulti)," genes with multiple unique integrations"))

  spIntronAll = split(usable,usable$MANE_intron)
  spIntron = spIntronAll[sapply(spIntronAll,nrow)>1]
  print(paste(length(spIntron)," introns with multiple integrations"))

  spCompositePos = split(usable,usable$CompositePos)
  spCompositePos = spCompositePos[sapply(spCompositePos,nrow)>1]
  # Collapse near-duplicate barcodes (likely PCR/sequencing error variants of the same
  # true barcode, not independent integrations) within each position before comparing.
  spCompositePos = lapply(spCompositePos, dedupBarcodesHamming, minDist = 3)
  spCompositePos = spCompositePos[sapply(spCompositePos,nrow)>1]
  print(paste(length(spCompositePos)," positions with multiple integrations (after Hamming-dist<3 barcode dedup)"))

  # Barcodes-per-gene and barcodes-per-intron use the full grouping (every gene/intron,
  # including singletons). Barcodes-per-position is restricted to the shared/multi-
  # barcode position subset (spCompositePos, same one used for the "Same position"
  # distance comparison above) -- almost every position has exactly 1 barcode, so an
  # unrestricted histogram would just be a spike at 1 and hide the (real) multi-barcode
  # tail we're actually interested in here.
  # A group can legitimately have zero members (e.g. no shared positions survive Hamming
  # dedup) -- sapply() on an empty list returns list(), not an empty numeric vector, which
  # otherwise makes data.frame() choke trying to recycle it against the length-1 "category".
  countDf = function(lst, label) {
    counts = if (length(lst)==0) numeric(0) else sapply(lst, nrow)
    data.frame(count = counts, category = rep(label, length(counts)))
  }
  histDf = rbind(
    countDf(spGene, "Barcodes per gene"),
    countDf(spIntronAll, "Barcodes per intron"),
    countDf(spCompositePos, "Barcodes per position\n(shared positions only)")
  )
  histDf$category = factor(histDf$category, levels = unique(histDf$category))
  pHist = ggplot(histDf, aes(x = count)) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "steelblue", color = "white") +
    facet_wrap(~category, scales = "free", ncol = 3) +
    labs(title = mapName, x = "Number of barcodes", y = "Count") +
    theme_bw(base_size = 16)
  print(pHist)

  # Select introns with multiple unique integrations
  spIntronMulti=spIntron[sapply(spIntron,function(x){length(unique(x$Pos))>1 & length(unique(x$Pos))==nrow(x)})]
  print(paste(length(spIntronMulti)," introns with multiple unique integrations"))

  inGene = usable[!is.na(usable$MANE_all),]
  nRandomPairs = 10000
  # BUG (memory blowout, fixed): this used to be dist(inGene[sample(...,nRandomPairs),"av"])
  # -- sampling nRandomPairs=10000 *barcodes* and then taking dist() (every pairwise
  # distance among them) is C(10000,2) ~= 50 million pairs, not the ~10,000-pair "Random
  # pairs" the label/variable name describe. That 50M-element vector then fed
  # wilcox.test(), ggplot's boxDf, etc. -- present in the unmodified original too (confirmed
  # by reproducing the same OOM kill with the pre-strand-fix script against current data),
  # not something introduced by the CompositePos strand-key change above. Sample
  # nRandomPairs pairs of barcodes directly instead (with replacement) so this stays O(10^4),
  # matching every other category's scale.
  randomPairDists = abs(inGene$av[sample(1:nrow(inGene), nRandomPairs, replace = TRUE)] -
                         inGene$av[sample(1:nrow(inGene), nRandomPairs, replace = TRUE)])
  # "Same gene" pairs would otherwise include pairs that are *also* "Same intron" pairs
  # (confirmed: ~28% of within-gene pairs share an intron here) -- exclude those so the
  # two categories are disjoint; genes left with no different-intron pair drop out (NA).
  intraGeneDists = sapply(spGeneMulti, meanDistDiffIntron)
  intraGeneDists = intraGeneDists[!is.na(intraGeneDists)]
  intraIntronDists = if (length(spIntronMulti)==0) numeric(0) else sapply(spIntronMulti,function(x){mean(dist(x$av))})
  intraPosDists = if (length(spCompositePos)==0) numeric(0) else sapply(spCompositePos,function(x){mean(dist(x$av))})

  boxData = list("Random\npairs" = randomPairDists,
                 "Same gene,\ndifferent introns" = intraGeneDists,
                 "Same\nintron" = intraIntronDists,
                 "Same\nposition" = intraPosDists)
  # Number of individual barcodes underlying each category (not the number of pairwise
  # distances, and not the number of values plotted -- each Same-x value is already a
  # per-group *mean* distance, one per gene/intron/position). For "Same gene, different
  # introns", only count barcodes in genes that actually contributed a valid
  # different-intron pair (i.e. survived the intraGeneDists NA-filter above). Alongside
  # that, also show how many distinct genes/introns/positions those barcodes fall into --
  # always <= the barcode count, since each group needs >=2 barcodes to contribute a
  # distance at all. "Random pairs" has no such grouping (each value is already one pair,
  # not barcodes pooled by gene/intron/position), so it's counted in pairs, not barcodes.
  sumNrow = function(lst) if (length(lst)==0) 0 else sum(sapply(lst, nrow))
  nBarcodes = c(nRandomPairs,
               sumNrow(spGeneMulti[names(intraGeneDists)]),
               sumNrow(spIntronMulti),
               sumNrow(spCompositePos))
  nGroups = c(NA, length(intraGeneDists), length(spIntronMulti), length(spCompositePos))
  groupUnit = c(NA, "genes", "introns", "positions")
  boxUnit = c("pairs", "barcodes", "barcodes", "barcodes")
  boxLabels = paste0(names(boxData), "\nn=", nBarcodes, " ", boxUnit,
                     ifelse(is.na(nGroups), "", paste0("\n(", nGroups, " ", groupUnit, ")")))
  boxDf = do.call(rbind, lapply(names(boxData), function(k) if (length(boxData[[k]])==0) NULL else data.frame(group = k, value = boxData[[k]])))
  boxDf$group = factor(boxDf$group, levels = names(boxData))
  comparisonsList = list(c("Random\npairs", "Same gene,\ndifferent introns"),
                         c("Same gene,\ndifferent introns", "Same\nintron"),
                         c("Same\nintron", "Same\nposition"))
  boxColors = setNames(c("grey70", colorRampPalette(c("#B2EBF2", "#006064"))(3)), names(boxData))
  boxTitle = if (nzchar(titleExtra)) paste(mapName, titleExtra) else mapName
  boxYLab = if (nzchar(titleExtra)) paste("Average distance of", titleExtra) else "Average distance"
  p = ggplot(boxDf, aes(x = group, y = value, fill = group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 20, outlier.size = 0.8) +
    stat_compare_means(comparisons = comparisonsList, method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = boxColors, guide = "none") +
    scale_x_discrete(labels = boxLabels) +
    labs(title = boxTitle, x = NULL, y = boxYLab) +
    theme_bw(base_size = 16)
  print(p)

  # ---------------------------------------------------------------------------
  # NEW: exp.log10 *level* distributions (not pairwise distances) across the same four
  # groupings -- a random sample of individual barcodes, and every barcode pooled from
  # every gene/intron/(Hamming-deduped) position that has >=2 barcodes. Unlike the
  # distance boxplot above (one value per gene/intron/position -- how spread out a
  # group's members are), this compares the actual exp.log10 values between barcodes in
  # each category, so it answers a different question: do barcodes that happen to share
  # a gene/intron/position tend to be higher- or lower-expressed than a random barcode?
  # No "different introns" exclusion here (unlike the distance boxplot's disjoint
  # categories above) -- a barcode legitimately counts toward both "share a gene" and
  # "share an intron" if applicable, since these are independent value-distribution
  # comparisons, not a partition.
  # ---------------------------------------------------------------------------
  levelsData = list(
    "Random\n10,000 barcodes" = inGene$av[sample(1:nrow(inGene), min(nRandomPairs, nrow(inGene)))],
    "Share a\ngene" = unlist(lapply(spGeneMulti, function(x) x$av), use.names = FALSE),
    "Share an\nintron" = unlist(lapply(spIntronMulti, function(x) x$av), use.names = FALSE),
    "Share a\nposition" = unlist(lapply(spCompositePos, function(x) x$av), use.names = FALSE)
  )
  nBarcodesLevels = sapply(levelsData, length)
  levelsLabels = paste0(names(levelsData), "\nn=", nBarcodesLevels, " barcodes")
  levelsDf = do.call(rbind, lapply(names(levelsData), function(k) if (length(levelsData[[k]])==0) NULL else data.frame(group = k, value = levelsData[[k]])))
  levelsDf$group = factor(levelsDf$group, levels = names(levelsData))
  levelsComparisons = list(c("Random\n10,000 barcodes", "Share a\ngene"),
                            c("Share a\ngene", "Share an\nintron"),
                            c("Share an\nintron", "Share a\nposition"))
  levelsColors = setNames(c("grey70", colorRampPalette(c("#FFCCBC", "#BF360C"))(3)), names(levelsData))
  levelsTitle = if (nzchar(titleExtra)) paste(mapName, titleExtra) else mapName
  levelsYLab = if (nzchar(titleExtra)) paste("exp.log10 of", titleExtra) else "exp.log10"
  pLevels = ggplot(levelsDf, aes(x = group, y = value, fill = group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 20, outlier.size = 0.8) +
    stat_compare_means(comparisons = levelsComparisons, method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = levelsColors, guide = "none") +
    scale_x_discrete(labels = levelsLabels) +
    labs(title = levelsTitle, x = NULL, y = levelsYLab) +
    theme_bw(base_size = 16)
  print(pLevels)

  # ---------------------------------------------------------------------------
  # NEW: % of exp.log10 variance explained by gene/intron/position identity (one-way
  # ANOVA eta-squared = SS_between/SS_total), restricted to the same >=2-barcode groups
  # used above. Each ANOVA's own F-test p-value is the appropriate significance test for
  # whether that grouping explains a real (non-random) share of the variance.
  # ---------------------------------------------------------------------------
  veGene = varianceExplained(spGeneMulti)
  veIntron = varianceExplained(spIntronMulti)
  vePos = varianceExplained(spCompositePos)
  veDf = data.frame(
    category = factor(c("Gene", "Intron", "Position"), levels = c("Gene", "Intron", "Position")),
    pct = c(veGene["pct"], veIntron["pct"], vePos["pct"]),
    p = c(veGene["p"], veIntron["p"], vePos["p"]),
    n = c(veGene["n"], veIntron["n"], vePos["n"]),
    nGroups = c(veGene["nGroups"], veIntron["nGroups"], vePos["nGroups"])
  )
  veDf$label = ifelse(is.na(veDf$p), "NA",
                       paste0("p=", formatC(veDf$p, format = "e", digits = 1),
                              "\nn=", veDf$n, " (", veDf$nGroups, " groups)"))
  veTitle = if (nzchar(titleExtra)) paste(mapName, titleExtra) else mapName
  veYLab = if (nzchar(titleExtra)) paste("% variance in", titleExtra, "explained") else "% variance in exp.log10 explained"
  veMax = suppressWarnings(max(veDf$pct, na.rm = TRUE))
  if (!is.finite(veMax)) veMax = 1
  pVarExplained = ggplot(veDf, aes(x = category, y = pct, fill = category)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_text(aes(label = label), vjust = -0.15, size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = c(Gene = "#4A90D9", Intron = "#F5A623", Position = "#D0021B"), guide = "none") +
    scale_y_continuous(limits = c(0, veMax * 1.3)) +
    labs(title = veTitle, x = NULL, y = veYLab) +
    theme_bw(base_size = 16)
  print(pVarExplained)

  # UpSet plot: which barcodes are expressed (exp.log10>0) in each library
  # BUG (memory blowout, fixed): min_size=0 doesn't add a floor -- it renders every
  # *observed* set-combination (up to 2^length(expNames)), not a fixed/small number. With
  # 12 libraries and ~18,700 barcodes here, most of the 4096 possible combinations are
  # observed at least once, so ComplexUpset was building/rendering a plot with thousands of
  # intersection columns -- this, not the CompositePos strand-key change or the earlier
  # random-pairs dist() bug, is what was still exhausting memory after those fixes (root-
  # caused via stdout buffering hiding how far the script actually got before OOM-kill: a
  # single-map, single-fix-at-a-time bisection pinned it here). n_intersections caps the
  # plot to the top-N largest intersections (standard UpSet convention), independent of how
  # many libraries/combinations exist.
  upsetDf = usable
  for (nameExp in expNames) {
    upsetDf[[nameExp]] = upsetDf[[paste0(nameExp, ".exp.log10")]] > 0
  }
  upsetDf = upsetDf[Reduce(`|`, upsetDf[expNames]), ]
  print(upset(upsetDf, intersect = expNames, name = "library", n_intersections = 40) +
          ggtitle(paste0("Expressed barcodes (exp.log10>0) - ", mapName)))
}
dev.off()
