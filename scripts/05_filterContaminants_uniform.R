contFactor = 30

path = "/home/projects/ulitsky/shared/data/human/TRIP/LP_run1//"
fileListFile = "LP.txt"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path=myArgs[1]
  fileListFile=myArgs[2]
  # Optional 3rd arg so the unified pipeline's config (CONT_FACTOR) can drive this without
  # editing the script; omit to keep the default above.
  if (length(myArgs)>=3) contFactor = as.numeric(myArgs[3])
}

tabs = list()
beds = list()
fileList = read.table(paste0(path,fileListFile))
for (name in fileList$V1)
{
  print(name)
  tabs[[name]]=read.table(paste0(path,name,".filtered_barcodes_poses.txt"),sep="\t")
  beds[[name]]=read.table(paste0(path,name,".filtered_barcodes.bed"),sep="\t")
  if (ncol(tabs[[name]]) != 11) stop(paste(name, "has", ncol(tabs[[name]]), "columns, expected 11"))
}

allBarcodes = list()
for (group in unique(fileList$V2))
{
  # Find the maximum count for each barcode
  allBarcodes[[group]] = c()
  for (name in fileList[fileList$V2==group,"V1"])
  {
    print(name)
    agg=aggregate(tabs[[name]][,10],by=tabs[[name]][,c(1,2,9)],FUN=max)
    allBarcodes[[group]][paste0(agg$V1,":",agg$V2,":",agg$V9)]=
      apply(data.matrix(cbind(agg$x,allBarcodes[[group]][paste0(agg$V1,":",agg$V2,":",agg$V9)])),1,max,na.rm=T)
  }
  print(paste(length(allBarcodes[[group]])," barcodes recoded"))
}
# Remove contaminants
for (name in fileList$V1)
{
  print(name)
  myGroup = fileList[fileList$V1==name,"V2"]
  allOtherBarcodes = c()
  for (otherGroup in setdiff(unique(fileList$V2),myGroup))
  {
    if (is.null(allOtherBarcodes[names(allBarcodes[[otherGroup]])]))
    {
      allOtherBarcodes=allBarcodes[[otherGroup]]
    } else {
      allOtherBarcodes[names(allBarcodes[[otherGroup]])]=
        apply(data.matrix(cbind(allBarcodes[[otherGroup]],allOtherBarcodes[names(allBarcodes[[otherGroup]])])),1,max,na.rm=T)
    }
  }
  # TRUE/FALSE/NA per (Chr,Pos,Barcode) key: keep if this key's own count is within
  # contFactor of the max count seen for the same key in any OTHER group; NA (kept) if the
  # key was never observed outside this group at all.
  keepIfEnough = function(keys, counts) {
    if (length(allOtherBarcodes)==0) return(rep(NA, length(keys)))
    counts >= (allOtherBarcodes[match(keys, names(allOtherBarcodes))]/contFactor)
  }

  myFilter = keepIfEnough(paste0(tabs[[name]]$V1,":",tabs[[name]]$V2,":",tabs[[name]]$V9), tabs[[name]][,10])
  print(paste(format(length(myFilter),big.mark=","),"barcodes considered"))
  print(paste(format(sum(is.na(myFilter)),big.mark=",")," not observed outside group"))
  print(paste(format(sum(myFilter,na.rm=T),big.mark=","),"to be retained"))
  print(paste(format(sum(!myFilter,na.rm=T),big.mark=","),"to be removed as contaminants"))
  head(myFilter[which(!myFilter)])
  write.table(subset(tabs[[name]],is.na(myFilter) | myFilter),paste0(path,name,".retained_barcodes_counts.txt"),col.names=F,row.names=F,quote=F,sep="\t")

  # .filtered_barcodes.bed is NOT row-aligned with .filtered_barcodes_poses.txt -- despite
  # often having the same row count, their rows differ in content (poses.txt applies
  # further position clustering/deduplication than the bed does), and a single bed row's
  # name field (V4) can itself list several near-duplicate barcode variants at one position
  # as "bc1:count1,bc2:count2,...". Reusing myFilter positionally (the original code) was
  # therefore filtering essentially unrelated rows -- confirmed on real data: a barcode
  # with 1 read here vs 755 in another group was correctly dropped from
  # .retained_barcodes_counts.txt but incorrectly survived in .retained_barcodes.bed.
  # Filter each bed row's sub-barcodes independently by their own (Chr,Start,Barcode) key,
  # keeping only the surviving ones (dropping the row entirely if none survive) rather than
  # assuming any correspondence with tabs[[name]]'s row order.
  bedRows = beds[[name]]
  splitParts = strsplit(as.character(bedRows$V4), ",", fixed=TRUE)
  partsPerRow = lengths(splitParts)
  flatParts = unlist(splitParts)
  rowIdx = rep(seq_len(nrow(bedRows)), partsPerRow)
  # Some rows' V4 is itself truncated mid-entry (a pre-existing upstream data issue, not
  # something introduced here -- confirmed on real data: a row whose V5 score claims 24
  # barcode variants only has ~10 complete ones before cutting off mid-barcode with no
  # ":count" suffix at all). A truncated fragment isn't a real barcode -- drop it outright
  # rather than defaulting it to "keep" the way a genuine "key not observed outside group"
  # NA is handled below.
  flatCount = suppressWarnings(as.numeric(sub("^.*:", "", flatParts)))
  malformed = is.na(flatCount)
  flatBarcode = sub(":[0-9]+$", "", flatParts)
  flatKeys = paste0(bedRows$V1[rowIdx], ":", bedRows$V2[rowIdx], ":", flatBarcode)
  flatKeep = keepIfEnough(flatKeys, flatCount)
  flatKeep[is.na(flatKeep) & !malformed] = TRUE
  flatKeep[malformed] = FALSE

  survivingParts = split(flatParts[flatKeep], rowIdx[flatKeep])
  survivingRows = as.integer(names(survivingParts))
  newV4 = character(nrow(bedRows))
  newV4[survivingRows] = sapply(survivingParts, paste, collapse=",")
  newScore = integer(nrow(bedRows))
  newScore[survivingRows] = lengths(survivingParts)
  bedRows$V4 = newV4
  bedRows$V5 = newScore
  bedFiltered = subset(bedRows, V4 != "")

  print(paste(format(sum(malformed),big.mark=","),"truncated/malformed bed barcode entries dropped (pre-existing upstream data issue)"))
  print(paste(format(sum(!flatKeep & !malformed),big.mark=","),"individual bed barcode entries removed as contaminants"))
  print(paste(format(nrow(bedRows)-nrow(bedFiltered),big.mark=","),"bed rows dropped entirely (no surviving barcode)"))
  write.table(bedFiltered,paste0(path,name,".retained_barcodes.bed"),col.names=F,row.names=F,quote=F,sep="\t")
}

