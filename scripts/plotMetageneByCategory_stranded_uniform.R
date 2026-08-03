#####
# Sibling of plotMetageneByCategory_singleCombined.R for strand-specific bigWig pairs (e.g.
# PRO-seq): takes a +strand and a -strand bigWig (same prefix, "<prefix>.plus.bw" /
# "<prefix>.minus.bw", e.g. PRO_Lis -> PRO_Lis.plus.bw/PRO_Lis.minus.bw) instead of one
# strand-agnostic bigWig. Barcodes on the "+" strand pull their signal from the plus
# bigWig, "-" strand barcodes from the minus bigWig (each strand's *own* track, not the
# genomic +/- of a single unstranded file) -- minus-strand bigWigs conventionally store
# values as negative (confirmed on PRO_Lis.minus.bw: range -15.5..0) purely so +/- tracks
# can be overlaid with opposite sign in a genome browser, so those values are negated back
# to a plain magnitude here before anything is compared/averaged. Both strands are then
# unified into one signal matrix at the metagene stage -- mean per position (not median: for
# a sparse point-signal assay like PRO-seq, checked empirically on real data that only ~53%
# of barcodes have *any* nonzero signal anywhere in their whole +-500bp window, so the
# per-position median -- and even the 75th percentile -- comes out exactly 0 at every single
# position; median/IQR is only meaningful for broad-coverage tracks like the non-stranded
# sibling's ATAC/H3K/POLR2A, not this kind of sparse data) -- see that script for the full
# annotation/categorization pipeline, unchanged here, including the 4 expression-level
# sub-lines (1<=exp.log10<2, 2<=exp.log10<3, 3<=exp.log10<4, exp.log10>=4) drawn within
# each of the four exp.log10>=1-gated category panels.
#####
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined/"
bwPath = "/home/projects/ulitsky/igoru/data/human/ENCODE/K562/"
mapNames = c("5ss", "CPA", "5ss_noP", "CPA_noP")
expListFile = "exp_list.txt"
bwPrefix = "PRO_Lis"
window = 500
smoothWindow = 10

myArgs = commandArgs(trailingOnly = TRUE)
if (length(myArgs) >= 1 && nzchar(myArgs[1])) mapNames = strsplit(myArgs[1], ",")[[1]]
if (length(myArgs) >= 2 && nzchar(myArgs[2])) expListFile = myArgs[2]
if (length(myArgs) >= 3 && nzchar(myArgs[3])) bwPrefix = myArgs[3]
print(mapNames)
print(expListFile)
print(bwPrefix)

catLevels = c("MANE_expressed", "GENCODE", "Antisense_expressed", "Intergenic_expressed", "Not expressed", "Other")
panelLevels = c("All", catLevels)
panelColors = setNames(c("black", brewer.pal(5, "Set1"), "gray50"), panelLevels)

# Categories whose defining rule already requires exp.log10>=1 -- these get split into 4
# expression-level lines per panel instead of one aggregate line.
expressedCats = c("MANE_expressed", "GENCODE", "Antisense_expressed", "Intergenic_expressed")
expBinLabels = c("1<=exp.log10<2", "2<=exp.log10<3", "3<=exp.log10<4", "exp.log10>=4")
expBinBreaks = c(1, 2, 3, 4, Inf)

# Centered moving average; leaves smoothWindow/2 NAs at each edge (a small, standard
# smoothing edge effect -- plot() just leaves a short gap there).
smoothSignal = function(x, window = smoothWindow) as.numeric(stats::filter(x, rep(1 / window, window), sides = 2))

bwPlus = BigWigFile(paste0(bwPath, bwPrefix, ".plus.bw"))
bwMinus = BigWigFile(paste0(bwPath, bwPrefix, ".minus.bw"))
chromSizes = seqlengths(bwPlus)  # confirmed identical to bwMinus's for PRO_Lis
relPos = -window:window

# ---- per-map: annotate, categorize, and average the (strand-flipped, strand-sourced)
# bigWig signal per panel -- returns a named list (by panelLevels) of
# list(n=, median=, q25=, q75=) for one map ------------------------------------------------
computeMapPanels = function(mapName)
{
  cat("\n====", mapName, "====\n")

  # quote="" and no explicit header= (autodetected via the ragged header -> row.names, same
  # as annotateTRIPGeneBodies_genuine_barcodes.R's read of *.annotated_tables.txt).
  annotTable = read.table(paste0(path, mapName, ".txt.annotated_table_used.txt"), sep = "\t", quote = "")
  annotTable$Barcode = rownames(annotTable)

  expList = read.table(paste0(path, expListFile), sep = "\t", header = FALSE,
                        col.names = c("mapType", "exp"))
  relevantExp = expList$exp[expList$mapType == mapName]
  if (length(relevantExp) == 0) stop(paste0("No experiments found for map '", mapName, "' in ", expListFile))
  cat("Experiments for", mapName, ":", paste(relevantExp, collapse = ", "), "\n")

  sumExp = NULL
  for (s in relevantExp) {
    tab = read.table(paste0(path, s, ".combined_counts.txt"), sep = "\t", header = TRUE)
    if (is.null(sumExp)) {
      sumExp = tab["exp"]
    } else {
      m = merge(sumExp, tab["exp"], by = 0, all = TRUE)
      rownames(m) = m$Row.names
      m[is.na(m)] = 0
      sumExp = data.frame(exp = m$exp.x + m$exp.y, row.names = rownames(m))
    }
  }
  sumExp$Barcode = rownames(sumExp)
  sumExp$exp.log10 = log10(1 + sumExp$exp)

  merged = merge(annotTable, sumExp[, c("Barcode", "exp", "exp.log10")], by = "Barcode")
  rownames(merged) = merged$Barcode

  # ---- mutually exclusive categories (first match wins; "Other" catches the rest) -------
  category = rep(NA_character_, nrow(merged))
  assignCat = function(cond, label) {
    hit = which(is.na(category) & cond)
    category[hit] <<- label
  }
  assignCat(!is.na(merged$MANE_all) & merged$MANE_all_exp.norm.log10 >= 1 & merged$exp.log10 >= 1, "MANE_expressed")
  assignCat(!is.na(merged$gencode_noncoding) & merged$exp.log10 >= 1, "GENCODE")
  assignCat(!is.na(merged$MANE_antisense) & merged$exp.log10 >= 1, "Antisense_expressed")
  assignCat(merged$exp.log10 >= 1, "Intergenic_expressed")
  assignCat(merged$exp.log10 <= 0.1, "Not expressed")
  assignCat(rep(TRUE, nrow(merged)), "Other")
  merged$Category = category

  cat("Barcode counts per category (", mapName, "):\n", sep = "")
  print(table(factor(merged$Category, levels = catLevels)))

  # ---- bigwig signal -500..+500 around each barcode: + strand barcodes from bwPlus,
  # - strand barcodes from bwMinus (sign-corrected), both windows already strand-flipped to
  # a consistent 5'->3' orientation -- then unified into one signal matrix -----------------
  used = merged
  gr = GRanges(seqnames = used$Chr,
               ranges = IRanges(start = used$Pos - window, end = used$Pos + window),
               strand = used$Strand)

  # import() errors (rather than silently clipping) on a range running off a chromosome end
  # or onto a seqname the bigwig doesn't have -- drop those few barcodes instead.
  inBounds = as.character(seqnames(gr)) %in% names(chromSizes) &
    start(gr) >= 1 & end(gr) <= chromSizes[as.character(seqnames(gr))]
  inBounds[is.na(inBounds)] = FALSE
  if (any(!inBounds)) cat(sum(!inBounds), "barcode window(s) dropped (off a chromosome end / unknown seqname)\n")
  gr = gr[inBounds]
  used = used[inBounds, ]

  sigMat = matrix(NA_real_, nrow = nrow(used), ncol = 2 * window + 1)

  plusIdx = which(as.character(strand(gr)) == "+")
  if (length(plusIdx) > 0) {
    sigListPlus = import(bwPlus, which = gr[plusIdx], as = "NumericList")
    for (i in seq_along(plusIdx)) {
      v = as.numeric(sigListPlus[[i]])
      v[is.na(v)] = 0
      sigMat[plusIdx[i], ] = v
    }
  }

  minusIdx = which(as.character(strand(gr)) == "-")
  if (length(minusIdx) > 0) {
    sigListMinus = import(bwMinus, which = gr[minusIdx], as = "NumericList")
    for (i in seq_along(minusIdx)) {
      v = as.numeric(sigListMinus[[i]])
      v[is.na(v)] = 0
      v = rev(-v)  # sign-correct (minus.bw stores negative magnitudes) and reverse orientation
      sigMat[minusIdx[i], ] = v
    }
  }

  # Per-position mean across barcodes, smoothed -- mean, not median (see header comment: for
  # a sparse assay like PRO-seq, the per-position median is 0 almost everywhere).
  meanSignal = function(idx) {
    if (length(idx) == 0) return(NULL)
    smoothSignal(colMeans(sigMat[idx, , drop = FALSE], na.rm = TRUE))
  }

  panels = lapply(panelLevels, function(p) {
    idx = if (p == "All") seq_len(nrow(used)) else which(used$Category == p)
    if (p %in% expressedCats) {
      expBin = cut(used$exp.log10[idx], breaks = expBinBreaks, labels = expBinLabels,
                   right = FALSE, include.lowest = TRUE)
      bins = lapply(expBinLabels, function(bl) {
        sub = idx[which(expBin == bl)]
        list(n = length(sub), signal = meanSignal(sub))
      })
      names(bins) = expBinLabels
      return(list(n = length(idx), signal = NULL, bins = bins))
    }
    list(n = length(idx), signal = meanSignal(idx), bins = NULL)
  })
  names(panels) = panelLevels
  panels
}

allPanels = lapply(mapNames, computeMapPanels)
names(allPanels) = mapNames

# ---- one combined PDF: rows = maps, columns = panels (All + each category) -------------
outFile = paste0(path, "plots/", paste(mapNames, collapse = "_"), ".", bwPrefix, ".metagene_by_category_grid.pdf")
pdf(outFile, width = 3.2 * length(panelLevels), height = 3.2 * length(mapNames), useDingbats = FALSE)
par(mfrow = c(length(mapNames), length(panelLevels)))
for (mapName in mapNames) {
  mapPanels = allPanels[[mapName]]
  sharedYlim = range(unlist(lapply(mapPanels, function(pd) {
    if (!is.null(pd$bins)) unlist(lapply(pd$bins, "[[", "signal")) else pd$signal
  })), na.rm = TRUE)
  for (p in panelLevels) {
    pd = mapPanels[[p]]
    if (!is.null(pd$bins)) {
      if (all(sapply(pd$bins, function(b) is.null(b$signal)))) {
        plot.new()
        title(main = paste0(mapName, " / ", p, " (n=0)"))
        next
      }
      shades = setNames(colorRampPalette(c("gray80", panelColors[p]))(4), expBinLabels)
      plot(relPos, rep(NA_real_, length(relPos)), type = "n", ylim = sharedYlim,
           main = paste0(mapName, " / ", p, " (n=", pd$n, ")"),
           xlab = "Position relative to barcode (nt)", ylab = paste(bwPrefix, "signal"))
      for (bl in expBinLabels) {
        if (is.null(pd$bins[[bl]]$signal)) next
        lines(relPos, pd$bins[[bl]]$signal, lwd = 2, col = shades[bl])
      }
      abline(v = 0, lty = 2, col = "gray50")
      legend("topright", legend = paste0(expBinLabels, " (n=", sapply(pd$bins, "[[", "n"), ")"),
             col = shades, lwd = 2, cex = 0.55, bty = "n")
    } else if (is.null(pd$signal)) {
      plot.new()
      title(main = paste0(mapName, " / ", p, " (n=0)"))
    } else {
      plot(relPos, pd$signal, type = "l", lwd = 2, col = panelColors[p], ylim = sharedYlim,
           main = paste0(mapName, " / ", p, " (n=", pd$n, ")"),
           xlab = "Position relative to barcode (nt)", ylab = paste(bwPrefix, "signal"))
      abline(v = 0, lty = 2, col = "gray50")
    }
  }
}
dev.off()
cat("wrote", outFile, "\n")
