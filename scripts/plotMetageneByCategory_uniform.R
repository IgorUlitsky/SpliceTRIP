#####
# For a comma-separated list of TRIP maps (XXX1,XXX2,... -- default the 4 single_combined
# maps: 5ss, CPA, 5ss_noP, CPA_noP), pools expression across the experiments listed for
# each map in exp_list.txt (2 tab-separated columns: mapType, exp sample name -- same
# sum(raw exp) across a map's listed experiments, then log10(1+exp), convention as
# plotMetageneCorrelations_singleCombined.R's getExpression(); merged onto
# <map>.txt.annotated_table_used.txt the same way annotateTRIPGeneBodies_genuine_barcodes.R
# merges combined_counts.txt onto the annotated table). Buckets barcodes into mutually
# exclusive expression/annotation categories (first-match-wins, in the order below;
# "Other" catches whatever matches none of the specific rules), prints the count in each,
# then plots a metagene (-500..+500 around each barcode, strand-flipped to a consistent
# 5'->3' orientation using the barcode's own Strand, 5% trimmed mean per position, smoothed
# with a 10nt moving average) of a K562 ENCODE bigwig track (YYY). Panels: "All" (every
# integration) plus each category; the four categories that require exp.log10>=1
# (MANE_expressed, GENCODE, Antisense_expressed, Intergenic_expressed) are further split
# into 3 expression-level lines within their own panel (1<=exp.log10<2, 2<=exp.log10<3,
# exp.log10>=3, light-to-dark shades of that panel's color), since "expressed" barcodes
# span a wide range of actual expression. One combined PDF: a grid
# with one row per map and one column per panel, y axis shared within each row (map).
#####
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined/"
bwPath = "/home/projects/ulitsky/igoru/data/human/ENCODE/K562/"
mapNames = c("5ss", "CPA", "5ss_noP", "CPA_noP")
expListFile = "exp_list.txt"
bwName = "ATAC"
window = 500
smoothWindow = 10
minBarcodesToPlot = 20  # panels/expression-bins with fewer barcodes than this are skipped

myArgs = commandArgs(trailingOnly = TRUE)
if (length(myArgs) >= 1 && nzchar(myArgs[1])) mapNames = strsplit(myArgs[1], ",")[[1]]
if (length(myArgs) >= 2 && nzchar(myArgs[2])) expListFile = myArgs[2]
if (length(myArgs) >= 3 && nzchar(myArgs[3])) bwName = myArgs[3]
print(mapNames)
print(expListFile)
print(bwName)

catLevels = c("MANE_expressed", "GENCODE", "Antisense_expressed", "Intergenic_expressed", "Not expressed", "Other")
panelLevels = c("All", catLevels)
panelColors = setNames(c("black", brewer.pal(5, "Set1"), "gray50"), panelLevels)

# Categories whose defining rule already requires exp.log10>=1 -- these get split into 4
# expression-level lines per panel instead of one aggregate line.
expressedCats = c("MANE_expressed", "GENCODE", "Antisense_expressed", "Intergenic_expressed")
expBinLabels = c("1<=exp.log10<2", "2<=exp.log10<3", "exp.log10>=3")
expBinBreaks = c(1, 2, 3, Inf)

# Centered moving average; leaves smoothWindow/2 NAs at each edge (a small, standard
# smoothing edge effect -- plot() just leaves a short gap there).
smoothSignal = function(x, window = smoothWindow) as.numeric(stats::filter(x, rep(1 / window, window), sides = 2))

bwFile = paste0(bwPath, bwName, ".K562.bw")
bw = BigWigFile(bwFile)
chromSizes = seqlengths(bw)
relPos = -window:window

# ---- per-map: annotate, categorize, and average the (strand-flipped) bigwig signal per
# panel -- returns a named list (by panelLevels) of either list(n=, signal=) (single line)
# or list(n=, bins=<one list(n=,signal=) per expBinLabels>) (expressedCats panels) --------
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

  # ---- bigwig signal -500..+500 around each barcode, strand-flipped to a consistent
  # 5'->3' orientation ---------------------------------------------------------------------
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

  sigList = import(bw, which = gr, as = "NumericList")
  sigMat = do.call(rbind, lapply(seq_along(sigList), function(i) {
    v = as.numeric(sigList[[i]])
    v[is.na(v)] = 0
    if (as.character(strand(gr)[i]) == "-") v = rev(v)
    v
  }))

  # 5% trimmed mean (drops the top/bottom 5% of barcodes at each position before
  # averaging) -- less sensitive to the single most extreme outlier barcode than a plain
  # mean, without collapsing to 0 the way the per-position median did (see header comment).
  # Below minBarcodesToPlot, skip the line entirely rather than plot it: with too few
  # barcodes, 5% trimming barely removes anything (n=20 -> 1 barcode off each end), so a
  # single outlier can still dominate the whole line.
  meanSignal = function(idx) {
    if (length(idx) < minBarcodesToPlot) return(NULL)
    smoothSignal(apply(sigMat[idx, , drop = FALSE], 2, mean, trim = 0.05, na.rm = TRUE))
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
outFile = paste0(path, "plots/", paste(mapNames, collapse = "_"), ".", bwName, ".metagene_by_category_grid.pdf")
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
      shades = setNames(colorRampPalette(c("gray80", panelColors[p]))(length(expBinLabels)), expBinLabels)
      plot(relPos, rep(NA_real_, length(relPos)), type = "n", ylim = sharedYlim,
           main = paste0(mapName, " / ", p, " (n=", pd$n, ")"),
           xlab = "Position relative to barcode (nt)", ylab = paste(bwName, "signal"))
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
           xlab = "Position relative to barcode (nt)", ylab = paste(bwName, "signal"))
      abline(v = 0, lty = 2, col = "gray50")
    }
  }
}
dev.off()
cat("wrote", outFile, "\n")
