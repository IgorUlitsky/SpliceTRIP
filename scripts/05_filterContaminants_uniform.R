# Stage 05: filterContaminants -- cross-group contamination removal, keyed PER BARCODE.
#
# THE BUG THIS FIXES
# ------------------
# The cross-group reference set (allBarcodes) was keyed on paste0(V1,":",V2,":",V9) and
# valued by V10. At a position where several near-duplicate barcode variants were collapsed
# together, V9 is not a single barcode -- it is the comma-joined LIST of every variant at
# that position, and V10 is the position-level MAX count, not any one barcode's count. Two
# groups almost never collapse the identical variant set (they differ with read depth), so
# the composite key strings differ, match() returns NA, keepIfEnough() returns NA, and the
# row is kept as "not observed outside group".
#
# Real example (run3 batch, contFactor 10), chr1:57942600:
#   5ss_Map_rep1_S81  V4 = TTTTCGAGATTGTACATC:2,CCTTCTGTTTGCTGGATC:1
#                     key = chr1:57942600:TTTTCGAGATTGTACATC,CCTTCTGTTTGCTGGATC
#   CPA_Map_rep2_S79  V4 = CCTTCTGTTTGCTGGATC:646,TTTTCGAGATTGTACATC:508,...(11 variants)
#                     key = chr1:57942600:CCTTCTGTTTGCTGGATC,TTTTCGAGATTGTACATC,...
# TTTTCGAGATTGTACATC has 2 reads in 5ss vs 508 in CPA -- a textbook contaminant (2 < 508/10)
# -- but the two composite keys never match, so it survived in both .retained_barcodes_counts.txt
# and .retained_barcodes.bed. ~2-12% of rows per library have a comma-joined V9, and those
# are exactly the high-depth positions where cross-library contamination concentrates.
#
# Note this also silently defeated the existing per-sub-barcode filtering of the .bed output:
# that code builds single-barcode keys (Chr:Start:barcode), which can never match a composite
# reference key, so every sub-barcode at a multi-variant position came back NA and was kept.
#
# THE FIX
# -------
# Build the reference set from V4 exploded into individual (barcode, count) pairs, keyed on
# Chr:Pos:barcode and valued by that barcode's OWN count. Both outputs are then filtered per
# sub-barcode against that reference. The .bed logic is unchanged -- it already worked this
# way and simply starts matching now. The counts/poses output is filtered per sub-barcode too
# (previously whole-row), rewriting V4/V7/V8/V9/V10/V11 to describe only the survivors.

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

# Explode a "bc1:count1,bc2:count2,..." column into one entry per barcode, carrying the
# originating row index so survivors can be reassembled. Entries that don't end in ":<digits>"
# are truncated/malformed upstream records (a pre-existing data issue: a row whose score
# claims N variants can cut off mid-barcode with no ":count" suffix) -- they are not real
# barcodes, so they are reported and dropped rather than defaulted to "keep".
explodeBarcodes = function(chrom, pos, field)
{
  parts = strsplit(as.character(field), ",", fixed=TRUE)
  rowIdx = rep(seq_along(field), lengths(parts))
  flat = unlist(parts)
  count = suppressWarnings(as.numeric(sub("^.*:", "", flat)))
  barcode = sub(":[0-9]+$", "", flat)
  list(rowIdx=rowIdx, part=flat, barcode=barcode, count=count,
       key=paste0(chrom[rowIdx], ":", pos[rowIdx], ":", barcode),
       malformed=is.na(count))
}

# Per-group reference: for every (Chr,Pos,barcode), the highest count that barcode reaches
# anywhere in the group. Keyed per barcode, NOT per collapsed variant-list -- this is the fix.
allBarcodes = list()
for (group in unique(fileList$V2))
{
  allBarcodes[[group]] = c()
  for (name in fileList[fileList$V2==group,"V1"])
  {
    print(name)
    e = explodeBarcodes(tabs[[name]]$V1, tabs[[name]]$V2, tabs[[name]]$V4)
    keep = !e$malformed
    if (!any(keep)) next
    # max count per key within this sample, then max against what the group already has.
    # tapply returns a named array -- rebuild as a plain named numeric so [<- extends it by
    # name. The empty-group case must be handled explicitly: indexing a NULL vector yields
    # NULL, and pmax(x, <zero-length>) is zero-length, which makes the assignment a silent
    # no-op (this left every group's reference set empty and disabled filtering entirely).
    perKeyArr = tapply(e$count[keep], e$key[keep], max)
    perKey = as.numeric(perKeyArr); names(perKey) = names(perKeyArr)
    if (length(allBarcodes[[group]]) == 0) {
      allBarcodes[[group]] = perKey
    } else {
      prev = allBarcodes[[group]][names(perKey)]
      prev[is.na(prev)] = -Inf
      allBarcodes[[group]][names(perKey)] = pmax(perKey, prev)
    }
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
  # TRUE/FALSE/NA per (Chr,Pos,Barcode) key: keep if this barcode's own count is within
  # contFactor of the max count that same barcode reaches at the same position in any OTHER
  # group; NA (kept) if the barcode was never observed at that position outside this group.
  keepIfEnough = function(keys, counts) {
    if (length(allOtherBarcodes)==0) return(rep(NA, length(keys)))
    counts >= (allOtherBarcodes[match(keys, names(allOtherBarcodes))]/contFactor)
  }

  # ---- counts/poses output, filtered per sub-barcode -------------------------------------
  tabRows = tabs[[name]]
  e = explodeBarcodes(tabRows$V1, tabRows$V2, tabRows$V4)
  eKeep = keepIfEnough(e$key, e$count)
  eKeep[is.na(eKeep) & !e$malformed] = TRUE
  eKeep[e$malformed] = FALSE

  print(paste(format(length(eKeep),big.mark=","),"barcode entries considered"))
  print(paste(format(sum(e$malformed),big.mark=","),"truncated/malformed entries dropped"))
  print(paste(format(sum(eKeep),big.mark=","),"to be retained"))
  print(paste(format(sum(!eKeep & !e$malformed),big.mark=","),"to be removed as contaminants"))

  # Reassemble surviving entries per row. V4/V9 keep their original order; V10/V11 are
  # recomputed by explicit max rather than taken from the first survivor -- V4 is usually but
  # NOT always sorted descending by count (confirmed on real data), so "first" is unsafe.
  survIdx = split(which(eKeep), e$rowIdx[eKeep])
  survRows = as.integer(names(survIdx))
  newV4 = character(nrow(tabRows)); newV9 = character(nrow(tabRows))
  newN  = integer(nrow(tabRows));   newV10 = tabRows$V10; newV11 = as.character(tabRows$V11)
  for (j in seq_along(survRows)) {
    ii = survIdx[[j]]; r = survRows[j]
    newV4[r] = paste(e$part[ii], collapse=",")
    newV9[r] = paste(e$barcode[ii], collapse=",")
    newN[r]  = length(ii)
    top = ii[which.max(e$count[ii])]
    newV10[r] = e$count[top]; newV11[r] = e$barcode[top]
  }
  tabRows$V4 = newV4; tabRows$V9 = newV9
  tabRows$V7 = newN; tabRows$V8 = newN      # V7/V8 are the variant count (verified on real data)
  tabRows$V10 = newV10; tabRows$V11 = newV11
  tabFiltered = subset(tabRows, V4 != "")
  print(paste(format(nrow(tabRows)-nrow(tabFiltered),big.mark=","),"poses rows dropped entirely (no surviving barcode)"))
  write.table(tabFiltered,paste0(path,name,".retained_barcodes_counts.txt"),col.names=F,row.names=F,quote=F,sep="\t")

  # ---- bed output ------------------------------------------------------------------------
  # .filtered_barcodes.bed is NOT row-aligned with .filtered_barcodes_poses.txt -- despite
  # often having the same row count, their rows differ in content (poses.txt applies
  # further position clustering/deduplication than the bed does), and a single bed row's
  # name field (V4) can itself list several near-duplicate barcode variants at one position
  # as "bc1:count1,bc2:count2,...". So each bed row's sub-barcodes are filtered independently
  # by their own (Chr,Start,Barcode) key rather than by any assumed row correspondence.
  bedRows = beds[[name]]
  be = explodeBarcodes(bedRows$V1, bedRows$V2, bedRows$V4)
  bKeep = keepIfEnough(be$key, be$count)
  bKeep[is.na(bKeep) & !be$malformed] = TRUE
  bKeep[be$malformed] = FALSE

  survingParts = split(be$part[bKeep], be$rowIdx[bKeep])
  survivingRows = as.integer(names(survingParts))
  newV4 = character(nrow(bedRows))
  newV4[survivingRows] = sapply(survingParts, paste, collapse=",")
  newScore = integer(nrow(bedRows))
  newScore[survivingRows] = lengths(survingParts)
  bedRows$V4 = newV4
  bedRows$V5 = newScore
  bedFiltered = subset(bedRows, V4 != "")

  print(paste(format(sum(be$malformed),big.mark=","),"truncated/malformed bed barcode entries dropped (pre-existing upstream data issue)"))
  print(paste(format(sum(!bKeep & !be$malformed),big.mark=","),"individual bed barcode entries removed as contaminants"))
  print(paste(format(nrow(bedRows)-nrow(bedFiltered),big.mark=","),"bed rows dropped entirely (no surviving barcode)"))
  write.table(bedFiltered,paste0(path,name,".retained_barcodes.bed"),col.names=F,row.names=F,quote=F,sep="\t")
}
