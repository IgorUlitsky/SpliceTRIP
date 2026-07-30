# =============================================================================
# Stage 18 — plotMetageneCorrelations
# =============================================================================
# Verbatim port of igoru's ~/R/trip/plotMetageneCorrelations_singleCombined.R, no logic
# changes (only the hardcoded default path/mapTypesList/expFile/outName at the top are
# this pipeline's usual per-script placeholders, overridden by args in practice). Reads
# stage 17's meta/<mapName>.combined_{histone,speckle}.txt matrices and correlates each
# barcode's summed expression (log10(1+exp), no norm) against every histone-mark and
# speckle/TSA-seq track's metagene signal -- per-construct bar charts plus a combined,
# p-value-annotated correlation heatmap across all constructs.
#
# Usage
#   Rscript 18_plotMetageneCorrelations_uniform.R <path> <corType> <mapTypesList> <expFile> <outName>
#
# Real example (as actually run for single_combined):
#   Rscript ~/R/trip/plotMetageneCorrelations_singleCombined.R \
#     /home/projects/ulitsky/shared/data/human/TRIP/single_combined/ pearson \
#     5ss,CPA,5ss_noP,CPA_noP exp_list.txt pearson_expression_singleCombined
#
# Inputs   histone_marks.txt, speckle_names.txt        (real, curated per-dataset files)
#          meta/<mapName>.combined_histone.txt          (stage 17)
#          meta/<mapName>.<firstHistoneMark>.tab.gz      (stage 17, read only for its
#                                                         barcode-name column)
#          meta/<mapName>.combined_speckle.txt          (stage 17)
#          meta/<mapName>.<firstSpeckleTrack>.tab.gz     (stage 17)
#          <expFile>                                    mapType<TAB>expName, one row per
#                                                        expName to sum for that construct
#                                                        (no norm -- replaces ### COMBOS
#                                                        entirely for this stage)
#          <expName>.combined_counts.txt                 (stage 07)
# Outputs  plots/<outName>.meta_corrs.pdf                per-construct histone/speckle
#                                                        correlation bar charts
#          plots/<outName>.summary_meta_corrs.pdf        combined heatmap (all constructs x
#                                                        all tracks), p<0.05 annotated
# =============================================================================

library(ComplexHeatmap)
library(circlize)
library(grid)

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined/"
corType = "pearson"
mapTypesList = "5ss,CPA,5ss_noP,CPA_noP"
expFile = "exp_list.txt"
outName = "pearson_expression_singleCombined"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs) > 0) {
  path = myArgs[1]
  corType = myArgs[2]
  mapTypesList = myArgs[3]
  expFile = myArgs[4]
  outName = myArgs[5]
}

myNames = strsplit(mapTypesList, ",")[[1]]

histoneNames = sort(read.table(paste0(path, "histone_marks.txt"))$V1)
speckleNames = sort(read.table(paste0(path, "speckle_names.txt"))$V1)

# expFile: 2 tab-separated columns (mapType, exp), one row per exp sample to include for
# that map type (multiple rows per map type are summed) -- replaces combos.txt entirely,
# no norm needed
expList = read.table(paste0(path, expFile), sep = "\t", header = F,
                      col.names = c("mapType", "exp"))

# ---- expression per barcode, per map type: sum(exp) only, log10(1+exp), no norm -------
getExpression = function(mapType) {
  expSamples = expList$exp[expList$mapType == mapType]
  stopifnot(length(expSamples) > 0)
  sumExp = NULL
  for (s in expSamples) {
    tab = read.table(paste0(path, s, ".combined_counts.txt"), sep = "\t", header = T)
    if (is.null(sumExp)) {
      sumExp = tab["exp"]
    } else {
      merged = merge(sumExp, tab["exp"], by = 0, all = T)
      rownames(merged) = merged$Row.names
      merged[is.na(merged)] = 0
      sumExp = data.frame(exp = merged$exp.x + merged$exp.y, row.names = rownames(merged))
    }
  }
  data.frame(exp.log10 = log10(1 + sumExp$exp), row.names = rownames(sumExp))
}

dfHistones = data.frame(row.names = histoneNames)
dfSpeckles = data.frame(row.names = speckleNames)
dfHistonesP = data.frame(row.names = histoneNames)
dfSpecklesP = data.frame(row.names = speckleNames)

pdf(paste0(path, "plots/", outName, ".meta_corrs.pdf"), width = 15)

for (name in myNames) {
  print(name)

  expTab = getExpression(name)

  combinedTab = read.table(paste0(path, "meta/", name, ".combined_histone.txt"), sep = "\t")
  names(combinedTab) = histoneNames
  sampleTab = read.table(gzfile(paste0(path, "meta/", name, ".", histoneNames[1], ".tab.gz")), sep = "\t")
  rownames(combinedTab) = sampleTab$V1

  combinedSpeckleTab = read.table(paste0(path, "meta/", name, ".combined_speckle.txt"), sep = "\t")
  names(combinedSpeckleTab) = speckleNames
  sampleHG19Tab = read.table(gzfile(paste0(path, "meta/", name, ".", speckleNames[1], ".tab.gz")), sep = "\t")
  rownames(combinedSpeckleTab) = sampleHG19Tab$V1

  merged = merge(expTab, combinedTab, by = 0)
  merged = merge(merged, combinedSpeckleTab, by.x = 1, by.y = 0)
  print(dim(merged))

  par(mar = c(10, 4, 4, 4))
  testsHistone = apply(merged[, histoneNames], 2, cor.test, merged$exp.log10, method = corType)
  corsHistone = sapply(testsHistone, "[[", "estimate")
  dfHistones[[name]] = corsHistone
  dfHistonesP[[name]] = sapply(testsHistone, "[[", "p.value")
  vec = sort(corsHistone)
  names(vec) = gsub(".rho$", "", names(vec))
  barplot(vec, las = 2, main = name, ylab = paste0(corType, " correlation"))

  par(mar = c(4, 25, 4, 4))
  testsSpeckles = apply(merged[, speckleNames], 2, cor.test, merged$exp.log10, method = corType)
  corsSpeckles = sapply(testsSpeckles, "[[", "estimate")
  dfSpeckles[[name]] = corsSpeckles
  dfSpecklesP[[name]] = sapply(testsSpeckles, "[[", "p.value")
  vec = sort(corsSpeckles)
  names(vec) = gsub(".rho$", "", names(vec))
  barplot(vec, las = 2, horiz = T, main = name, xlab = paste0(corType, " correlation"))
}
dev.off()

corLabel = paste0(toupper(substring(corType, 1, 1)), substring(corType, 2), " correlation")
fontsizeLabels = 8
legendParams = list(title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 14))

# page size sized to content so row/column labels never run off the page. Row-label width
# is measured directly via grid (not estimated from char count) since actual rendered text
# width is what matters; legend goes on the bottom (not right) so it never competes with
# row labels for horizontal space, which was overlapping them before.
maxTextWidthIn = function(labels, fontsize) {
  pdf(NULL)
  w = max(sapply(labels, function(l) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(l, gp = grid::gpar(fontsize = fontsize))), "inches", valueOnly = TRUE)
  }))
  dev.off()
  w
}
# NOTE: ComplexHeatmap's row_names_max_width defaults to a fixed unit(6, "cm") *regardless*
# of how much device space is available, so row names silently truncate at ~2.36in even on
# an oversized page unless this is explicitly overridden per-heatmap below.
estimatePdfSize = function(df, fontsize = fontsizeLabels, cellWidthIn = 0.45, cellHeightIn = 0.32) {
  # small safety margin: text-width measurement device vs. final render device can differ
  # by a hair, which is enough to clip the last character of the very longest label
  rowLabelWidthIn = maxTextWidthIn(rownames(df), fontsize) + 0.15
  bodyWidthIn = ncol(df) * cellWidthIn
  bodyHeightIn = nrow(df) * cellHeightIn
  width = bodyWidthIn + rowLabelWidthIn + 1
  height = bodyHeightIn + 4
  list(width = width, height = height, rowLabelWidthIn = rowLabelWidthIn,
       bodyWidthIn = bodyWidthIn, bodyHeightIn = bodyHeightIn)
}
szH = estimatePdfSize(dfHistones)
szS = estimatePdfSize(dfSpeckles)
pdfWidth = max(szH$width, szS$width)
pdfHeight = max(szH$height, szS$height)

# show the p-value (sci. notation, %3.2g) in cells where it's significant (<0.05)
pvalCellFun = function(pmat) {
  function(j, i, x, y, width, height, fill) {
    p = pmat[i, j]
    if (!is.na(p) && p < 0.05) {
      grid.text(sprintf("%3.2g", p), x, y, gp = gpar(fontsize = 5))
    }
  }
}

# symmetric around 0 (white), scaled to the largest |correlation| actually observed
# (shared across both heatmaps so their color scales stay comparable to each other)
maxAbsCor = max(abs(c(as.matrix(dfHistones), as.matrix(dfSpeckles))), na.rm = TRUE)
corColors = colorRamp2(c(-maxAbsCor, 0, maxAbsCor), c("blue", "white", "red"))

pdf(paste0(path, "plots/", outName, ".summary_meta_corrs.pdf"), width = pdfWidth, height = pdfHeight)
draw(Heatmap(as.matrix(dfHistones), name = corLabel, col = corColors,
             width = unit(szH$bodyWidthIn, "in"), height = unit(szH$bodyHeightIn, "in"),
             row_names_max_width = unit(szH$rowLabelWidthIn, "in"),
             row_names_gp = gpar(fontsize = fontsizeLabels), column_names_gp = gpar(fontsize = fontsizeLabels),
             heatmap_legend_param = legendParams,
             cell_fun = pvalCellFun(as.matrix(dfHistonesP))),
     heatmap_legend_side = "bottom")
draw(Heatmap(as.matrix(dfSpeckles), name = corLabel, col = corColors,
             width = unit(szS$bodyWidthIn, "in"), height = unit(szS$bodyHeightIn, "in"),
             row_names_max_width = unit(szS$rowLabelWidthIn, "in"),
             row_names_gp = gpar(fontsize = fontsizeLabels), column_names_gp = gpar(fontsize = fontsizeLabels),
             heatmap_legend_param = legendParams,
             cell_fun = pvalCellFun(as.matrix(dfSpecklesP))),
     heatmap_legend_side = "bottom")
dev.off()
