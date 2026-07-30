# =============================================================================
# Stage 16 — compareReplicatesAndMultiIntegrations
# =============================================================================
# Verbatim port of igoru's ~/R/trip/compareReplicatesAndMultiIntegrations_uniform.R, no
# logic changes. QC report per construct: how much do independent barcodes landing in the
# same gene / same intron / same exact position agree with each other, compared to random
# pairs? (expression distance boxplots with Wilcoxon comparisons), how many barcodes share
# a gene/intron/position (histograms), and which libraries express which barcodes (UpSet
# plot). Same reads-annotated_tables.txt-plus-combined_counts.txt-per-expName merge pattern
# as stage 11/12/13, but the actual comparisons (gene/intron/position multi-integration
# agreement) are unique to this stage.
#
# Usage
#   Rscript 16_compareReplicatesAndMultiIntegrations_uniform.R <path> <mapName1,mapName2,...> <expListFile> <outFile.pdf> [titleExtra]
#
# Real example (as actually run for IntronTRIP_combined):
#   bsub -R rusage[mem=30000] Rscript ~/R/trip/compareReplicatesAndMultiIntegrations_uniform.R \
#     ./ A1,B1,C1,D1,E1,B3-6,B3-9,H1 expression.exp_list.txt \
#     plots/replicates_integrations_IntronTRIP_expression.pdf
#
# expListFile: mapName<TAB>expName, one row per (construct, expression-readout) pair to
# include -- like stage 11, real production expListFiles are often hand-curated and NOT 1:1
# with ### COMBOS (e.g. only the "Total" readout per rep, as in the real
# expression.exp_list.txt this stage was designed around). See examples/exp_lists/ for real
# ones; pass your own expListFile as the 3rd arg to run this directly rather than through
# run_pipeline.sh if you need something ### COMBOS can't auto-generate.
#
# Inputs (per mapName)   <mapName>.txt.annotated_tables.txt              (stage 10)
# Inputs (per expName)   <expName>.combined_counts.txt                   (stage 07)
# Output                 one multi-page PDF (outFile), 3 pages per mapName: barcodes-per-
#                        gene/intron/position histograms, same-gene/intron/position vs
#                        random-pairs expression-distance boxplots, and an UpSet plot of
#                        which libraries express which barcodes.
# =============================================================================

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

  usable$CompositePos = paste0(usable$Chr,":",usable$Pos)

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
  distMat = dist(inGene[sample(1:nrow(inGene),min(nRandomPairs,nrow(inGene))),"av"])
  # "Same gene" pairs would otherwise include pairs that are *also* "Same intron" pairs
  # (confirmed: ~28% of within-gene pairs share an intron here) -- exclude those so the
  # two categories are disjoint; genes left with no different-intron pair drop out (NA).
  intraGeneDists = sapply(spGeneMulti, meanDistDiffIntron)
  intraGeneDists = intraGeneDists[!is.na(intraGeneDists)]
  intraIntronDists = if (length(spIntronMulti)==0) numeric(0) else sapply(spIntronMulti,function(x){mean(dist(x$av))})
  intraPosDists = if (length(spCompositePos)==0) numeric(0) else sapply(spCompositePos,function(x){mean(dist(x$av))})

  boxData = list("Random\npairs" = as.numeric(distMat),
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
  # distance at all. "Random pairs" has no such grouping (it's individual sampled
  # barcodes, not barcodes pooled by gene/intron/position), so it only gets the one count.
  sumNrow = function(lst) if (length(lst)==0) 0 else sum(sapply(lst, nrow))
  nBarcodes = c(min(nRandomPairs,nrow(inGene)),
               sumNrow(spGeneMulti[names(intraGeneDists)]),
               sumNrow(spIntronMulti),
               sumNrow(spCompositePos))
  nGroups = c(NA, length(intraGeneDists), length(spIntronMulti), length(spCompositePos))
  groupUnit = c(NA, "genes", "introns", "positions")
  boxLabels = paste0(names(boxData), "\nn=", nBarcodes, " barcodes",
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

  # UpSet plot: which barcodes are expressed (exp.log10>0) in each library
  upsetDf = usable
  for (nameExp in expNames) {
    upsetDf[[nameExp]] = upsetDf[[paste0(nameExp, ".exp.log10")]] > 0
  }
  upsetDf = upsetDf[Reduce(`|`, upsetDf[expNames]), ]
  print(upset(upsetDf, intersect = expNames, name = "library", min_size = 0) +
          ggtitle(paste0("Expressed barcodes (exp.log10>0) - ", mapName)))
}
dev.off()
