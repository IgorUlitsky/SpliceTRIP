# =============================================================================
# Stage 11 — pairwise shared-position comparison across mapping libraries
# =============================================================================
# Usage
#   Rscript 11_comparePairwiseSharedPositions_uniform.R <path> <mapName1,mapName2,...> <expListFile> <outFile.pdf> [titleExtra]
#
# Inputs   <mapName>.txt.annotated_tables.txt   (stage 10)
#          <expName>.combined_counts.txt        (stage 07, OR any other pipeline that
#                                                produces the same columns -- e.g. a
#                                                splicing-ratio exp.log10 instead of the
#                                                standard log10(1+exp), which this script
#                                                detects and preserves automatically; see
#                                                examples/exp_lists/ for real spliced/
#                                                unspliced-style expListFile examples)
#   expListFile: <mapName><TAB><expName>, one row per (construct, expression library) pair
#                to pool into that construct's per-position average -- see
#                examples/exp_lists/ for real ones, or auto-generate from ### COMBOS
#                (run_pipeline.sh's "annotated" target does this automatically).
#
# Outputs  <outFile.pdf>                                    4 ComplexHeatmaps (see below)
#          <outFile.pdf>.shared_positions.<mapA>_vs_<mapB>.txt   per unordered construct
#                                                            pair, one row per shared
#                                                            position
#
# For every pair of mapping libraries (every ### FASTQ_MAPS mapName, not every
# ### COMBOS row): how many genomic positions are shared (both libraries have >=1
# barcode mapped there), and further split into positions that also share an identical
# barcode (much stronger evidence of a true shared integration than merely landing on
# the same position independently) vs. ones that don't. For the shared positions, the
# Spearman correlation of the two libraries' per-position average exp.log10/map/norm,
# split the same way. Four ComplexHeatmaps: shared-position counts (top: total, bottom:
# with an identical barcode), then one each for exp.log10/map/norm correlation (top:
# different-barcode subset, bottom: same-barcode subset).
#
# No hardcoded paths in this script (unlike stages 08-10) -- everything it reads comes
# from other pipeline stages' own output, addressed relative to <path>.
# =============================================================================
library(ComplexHeatmap)
library(circlize)
library(grid)

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined///"

mapNames = c("5ss","CPA","5ss_noP","CPA_noP")
expListFile = "exp_list.txt"
outFile = "plots/pairwise_shared_positions_single.pdf"
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

# For each mapping library: load + merge exactly like
# compareReplicatesAndMultiIntegrations_uniform.R (same all.x=TRUE merge, same
# "exp.log10 already present" handling), then collapse from one row per barcode down to
# one value per genomic position -- for exp.log10, map, and norm alike -- by averaging
# each barcode's own per-expName-averaged value across every barcode mapped to that
# position within this library. Also keep the set of barcodes at each position (needed
# below for the "shared position with an identical barcode" split) and the total barcode
# count (shown in the heatmap axis labels).
posDataByMap = list()
for (mapName in mapNames)
{
  print(mapName)
  annotTab=read.table(paste0(path,mapName,".txt.annotated_tables.txt"),sep="\t",comment="")
  allMerged = annotTab
  expNames = expList[expList$V1==mapName,"V2"]
  for (nameExp in expNames)
  {
    countsTab =  read.table(paste0(path, nameExp, ".combined_counts.txt"), sep = "\t", quote = "")
    countsTab$Barcode = rownames(countsTab)
    annotTab$Barcode = rownames(annotTab)
    # If the file already has an exp.log10 column (e.g. a bounded splicing ratio), trust
    # it as-is rather than overwriting it with the generic log10(1+exp) transform.
    if (!("exp.log10" %in% names(countsTab))) {
      countsTab$exp.log10 = log(1+countsTab$exp,10)
    }
    countsTab$exp.norm = log((1 + countsTab$exp) / (1 + countsTab$norm), 2)
    for (col in c("norm","exp","map","exp.log10","exp.norm","exp_spliced","exp_total","exp_unspliced")) {
      if (col %in% names(countsTab)) names(countsTab)[names(countsTab)==col] = paste0(nameExp,".",col)
    }
    allMerged = merge(allMerged,countsTab,by=c("Barcode","Chr","Pos","Strand"),all.x=TRUE)
    allMerged$Row.names=NULL
  }
  allMerged$av=apply(allMerged[,grep("exp.log10$",names(allMerged)),drop=FALSE],1,mean,na.rm=T)
  # Named explicitly per expName (not a blanket ".map$"/".norm$" suffix regex) --
  # "{nameExp}.exp.norm" also ends in ".norm", so a suffix match would wrongly pull it in
  # alongside the raw "{nameExp}.norm" column.
  mapCols = intersect(paste0(expNames,".map"), names(allMerged))
  normCols = intersect(paste0(expNames,".norm"), names(allMerged))
  allMerged$av.map=apply(allMerged[,mapCols,drop=FALSE],1,mean,na.rm=T)
  allMerged$av.norm=apply(allMerged[,normCols,drop=FALSE],1,mean,na.rm=T)

  usable = subset(allMerged,!is.na(allMerged$av))
  usable$CompositePos = paste0(usable$Chr,":",usable$Pos,":",usable$Strand)

  aggMean = function(col) setNames(aggregate(usable[[col]], by=list(CompositePos=usable$CompositePos), FUN=mean, na.rm=TRUE)$x,
                                    aggregate(usable[[col]], by=list(CompositePos=usable$CompositePos), FUN=mean, na.rm=TRUE)$CompositePos)
  posDataByMap[[mapName]] = list(
    av = aggMean("av"),
    avMap = aggMean("av.map"),
    avNorm = aggMean("av.norm"),
    barcodes = split(usable$Barcode, usable$CompositePos),
    nBarcodes = nrow(usable)
  )
  print(paste(length(posDataByMap[[mapName]]$av)," positions in ",mapName))
}

# Axis labels showing each library's total (usable) barcode count, used on every heatmap.
axisLabels = paste0(mapNames, "\n(", sapply(mapNames, function(m) posDataByMap[[m]]$nBarcodes), " barcodes)")

# Pairwise, for every off-diagonal (i,j): total shared positions, shared positions where
# the two libraries also have an identical barcode present, and (for those shared
# positions) the Spearman correlation of the two libraries' per-position averages -- for
# exp.log10, map, and norm each, split into the subset of shared positions that DO vs
# DON'T have an identical barcode (a real biological distinction: a shared identical
# barcode at a shared position is much stronger evidence of a true shared/duplicated
# integration event than two libraries' barcodes merely landing on the same position
# independently). All matrices are naturally symmetric (intersect(A,B)==intersect(B,A)).
# The diagonal is left NA/blank rather than showing self-vs-self: a library's own total
# position count is typically far larger than any cross-library shared count, which would
# otherwise dominate the count heatmap's color scale and wash out the off-diagonal cells
# that are the actual point of "every pair of mapping libraries".
#
# For each pair (written once per unordered pair, i<j) a text file listing every shared
# position, its shared barcode(s) (if any), and both libraries' per-position average
# exp.log10/map/norm -- the exact numbers the correlations above are computed from.
n = length(mapNames)
minSharedForCorr = 3
countMat = matrix(NA_integer_, n, n, dimnames=list(mapNames,mapNames))
identicalBarcodeMat = matrix(NA_integer_, n, n, dimnames=list(mapNames,mapNames))
corrMatSame = matrix(NA_real_, n, n, dimnames=list(mapNames,mapNames))
corrMatDiff = matrix(NA_real_, n, n, dimnames=list(mapNames,mapNames))
mapCorrMatSame = matrix(NA_real_, n, n, dimnames=list(mapNames,mapNames))
mapCorrMatDiff = matrix(NA_real_, n, n, dimnames=list(mapNames,mapNames))
normCorrMatSame = matrix(NA_real_, n, n, dimnames=list(mapNames,mapNames))
normCorrMatDiff = matrix(NA_real_, n, n, dimnames=list(mapNames,mapNames))

corrIfEnough = function(x, y) if (length(x) >= minSharedForCorr) cor(x, y, method="spearman") else NA_real_

for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) next
    di = posDataByMap[[mapNames[i]]]
    dj = posDataByMap[[mapNames[j]]]
    shared = intersect(names(di$av), names(dj$av))
    hasSameBarcode = sapply(shared, function(pos) length(intersect(di$barcodes[[pos]], dj$barcodes[[pos]])) > 0)
    sharedSame = shared[hasSameBarcode]
    sharedDiff = shared[!hasSameBarcode]

    countMat[i,j] = length(shared)
    identicalBarcodeMat[i,j] = sum(hasSameBarcode, na.rm=T)

    corrMatSame[i,j] = corrIfEnough(di$av[sharedSame], dj$av[sharedSame])
    corrMatDiff[i,j] = corrIfEnough(di$av[sharedDiff], dj$av[sharedDiff])
    mapCorrMatSame[i,j] = corrIfEnough(di$avMap[sharedSame], dj$avMap[sharedSame])
    mapCorrMatDiff[i,j] = corrIfEnough(di$avMap[sharedDiff], dj$avMap[sharedDiff])
    normCorrMatSame[i,j] = corrIfEnough(di$avNorm[sharedSame], dj$avNorm[sharedSame])
    normCorrMatDiff[i,j] = corrIfEnough(di$avNorm[sharedDiff], dj$avNorm[sharedDiff])

    if (i < j && length(shared) > 0) {
      sharedBarcodesStr = sapply(shared, function(pos) paste(intersect(di$barcodes[[pos]], dj$barcodes[[pos]]), collapse=","))
      outDf = data.frame(CompositePos = shared, SharedBarcodes = sharedBarcodesStr, stringsAsFactors = FALSE)
      outDf[[paste0(mapNames[i],".exp.log10")]] = di$av[shared]
      outDf[[paste0(mapNames[j],".exp.log10")]] = dj$av[shared]
      outDf[[paste0(mapNames[i],".map")]] = di$avMap[shared]
      outDf[[paste0(mapNames[j],".map")]] = dj$avMap[shared]
      outDf[[paste0(mapNames[i],".norm")]] = di$avNorm[shared]
      outDf[[paste0(mapNames[j],".norm")]] = dj$avNorm[shared]
      write.table(outDf, paste0(path, outFile, ".shared_positions.", mapNames[i], "_vs_", mapNames[j], ".txt"),
                  sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
}

pdf(paste0(path,outFile), width=10, height=9)

# draw() already renders to the open device as a side effect -- wrapping it in print()
# would render it a second time, producing a duplicate page.
countTitle = if (nzchar(titleExtra)) paste("Shared positions -", titleExtra) else "Shared positions"
countMax = max(countMat, na.rm = TRUE)
draw(Heatmap(countMat, name="shared\npositions",
             col=colorRamp2(c(0, countMax), c("white","steelblue")),
             na_col = "grey90",
             cluster_rows=FALSE, cluster_columns=FALSE,
             row_labels=axisLabels, column_labels=axisLabels,
             column_title=paste0(countTitle, "\n(top: positions shared; bottom: of those, sharing an identical barcode)"),
             cell_fun=function(j,i,x,y,width,height,fill) {
               if (is.na(countMat[i,j])) return()
               # grid y increases upward, so "+ offset" is the top line and "- offset" the
               # bottom one.
               grid.text(sprintf("%d", countMat[i,j]), x, y + unit(2.2,"mm"), gp=gpar(fontsize=11))
               grid.text(sprintf("%d", identicalBarcodeMat[i,j]), x, y - unit(2.2,"mm"), gp=gpar(fontsize=9, col="grey30"))
             }))

# Each cell shows two Spearman rho values stacked: top = shared positions with different
# barcodes, bottom = shared positions that also share an identical barcode. Either half can
# be NA (blank) if fewer than minSharedForCorr positions fall in that subset -- e.g. two
# genuinely independent barcode pools will have ~0 identical-barcode positions, so the
# bottom number is expected to be blank for those pairs.
drawCorrHeatmap = function(matSame, matDiff, valueLabel) {
  title = if (nzchar(titleExtra)) paste0("Spearman correlation of shared-position ", valueLabel, " (", titleExtra, ")") else paste0("Spearman correlation of shared-position ", valueLabel)
  title = paste0(title, "\n(top: different barcode; bottom: same barcode)")
  draw(Heatmap(matDiff, name="Spearman\nrho\n(diff. bc)",
               col=colorRamp2(c(-1,0,1), c("#2166AC","white","#B2182B")),
               na_col = "grey90",
               cluster_rows=FALSE, cluster_columns=FALSE,
               row_labels=axisLabels, column_labels=axisLabels,
               column_title=title,
               cell_fun=function(j,i,x,y,width,height,fill) {
                 if (!is.na(matDiff[i,j])) grid.text(sprintf("%.2f", matDiff[i,j]), x, y + unit(2.2,"mm"), gp=gpar(fontsize=10))
                 if (!is.na(matSame[i,j])) grid.text(sprintf("%.2f", matSame[i,j]), x, y - unit(2.2,"mm"), gp=gpar(fontsize=10, col="grey30"))
               }))
}
drawCorrHeatmap(corrMatSame, corrMatDiff, "exp.log10")
drawCorrHeatmap(mapCorrMatSame, mapCorrMatDiff, "map")
drawCorrHeatmap(normCorrMatSame, normCorrMatDiff, "norm")

dev.off()
