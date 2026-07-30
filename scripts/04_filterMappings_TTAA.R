library(DescTools)

# Load parameters
path = "/home/labs/ulitsky/shared/data/human/TRIP/IntronTRIP_run3//"
#path = "/home/labs/ulitsky/lidiag/TRIP_data/In"
#name= "Map_5ss_noP_rep2_S96"
name = "IntronTRIP_Map_A1_Pool_1_2_S63"
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path=myArgs[1]
  name=myArgs[2]
}

# ============================================================================
# ===== NEW BLOCK (TTAA filter) — argument parsing + reference loading ========
# Usage:
#   default, NO filtering :  Rscript filterMappings_TTAA.R  <path>  <name>
#   with TTAA filtering   :  Rscript filterMappings_TTAA.R  <path>  <name>  ttaaFilter  <ttaa_bed_path>
# <ttaa_bed_path> = a 6-col BED like hg38_ttaa.bed (col2 = 0-based TTAA start).
# Without the literal flag "ttaaFilter" the script behaves EXACTLY like the original.
ttaaFilter  = FALSE
ttaaByChrom = NULL
if (length(myArgs) >= 3 && myArgs[3] == "ttaaFilter") {
  ttaaFilter = TRUE
  if (length(myArgs) < 4) stop("ttaaFilter requested but no TTAA bed path given (4th argument)")
  ttaaBedPath = myArgs[4]
  print(paste("TTAA filter ON; loading", ttaaBedPath))
  ttaaBed = read.table(ttaaBedPath, sep="\t", header=F,
                       colClasses=c("character","integer","NULL","NULL","NULL","NULL"))
  colnames(ttaaBed) = c("chrom","start")                       # start = 0-based TTAA start
  ttaaByChrom = lapply(split(ttaaBed$start, ttaaBed$chrom), sort)
  print(paste("Loaded", nrow(ttaaBed), "TTAA sites on", length(ttaaByChrom), "contigs"))
} else {
  print("TTAA filter OFF (default behaviour)")
}
# ===== end NEW BLOCK ========================================================
# ============================================================================

print("Loading datasets")
countTab <- read.table(paste0(path,name,".barcodes_poses.txt"),sep="\t",header=F)
bedTab <- read.table(paste0(path,name,".barcodes.bed"),sep="\t",header=F)
dim(countTab)
dim(bedTab)

countTab$MaxCount = as.numeric(sapply(strsplit(sapply(strsplit(as.character(countTab$V4),","),"[",1),":"),"[",2))
countTab$TopBarcode = sapply(strsplit(sapply(strsplit(as.character(countTab$V4),","),"[",1),":"),"[",1)
spCountTab = split(countTab,countTab$V1)
spBedTab = split(bedTab,bedTab$V1)

allOutTab = data.frame(matrix(nrow = 0,ncol=length(countTab)))
allOutBED = data.frame(matrix(nrow = 0,ncol=length(bedTab)))

print("Filtering")
minDist = 10
minFrac = 0.25
tooClose = 0
legit = 0
nonTTAA = 0                                            # ===== NEW: counter of positions dropped as non-TTAA
topDists = c()
addToBin <-function(bin,x)
{
  if (is.null(bin))
  {
    df=data.frame(matrix(nrow=0,ncol=length(countTab)))
    names(df)=names(countTab)
  } else {
    df=bin
  }
  df[nrow(df)+1,]=x
  df
}
for (chr in names(spCountTab))
{
  print(chr)
  chrTab = spCountTab[[chr]]
  chrBed = spBedTab[[chr]]
  chrTab=chrTab[order(chrTab$V2,chrTab$V6),]
  chrBed=chrBed[order(chrBed$V2,chrBed$V6),]

  # ===== NEW BLOCK (TTAA filter) — per-row TTAA membership on the SORTED chrTab
  # + strand: TTAA starts at V2-4 ;  - strand: TTAA starts at V2+1   (0-based, strand-aware).
  # ttaaOK is index-aligned to chrTab, i.e. to loop variable i below. Exact match, no window.
  # NOTE (minus-strand offset, V2+1 not V2): this pipeline's stage 03 records the
  # minus-strand position as V3-1 (the "V3-1" convention -- the last aligned base, since
  # BED end is exclusive), not V3. Empirically verified against the real hg38 TTAA
  # reference on real data: with a V3-1-produced barcodes.bed, tt=V2 gives a 0.09% hit
  # rate (i.e. broken) while tt=V2+1 gives ~94% -- matching the "+" strand's own hit rate
  # and confirming V2+1 recovers the correct site. If stage 03 is ever swapped back to the
  # uncorrected V3 convention, this must revert to tt=V2 (its own matching offset).
  if (ttaaFilter) {
    vec = ttaaByChrom[[chr]]                                    # sorted 0-based TTAA starts (NULL if contig absent)
    if (is.null(vec)) {
      print(paste("  WARNING: no TTAA sites for", chr, "-> all its positions will be dropped"))
      ttaaOK = rep(FALSE, nrow(chrTab))
    } else {
      tt = ifelse(chrTab$V6 == "+", chrTab$V2 - 4L, chrTab$V2 + 1L) # expected TTAA start
      ttaaOK = tt %in% vec
    }
  } else {
    ttaaOK = rep(TRUE, nrow(chrTab))
  }
  # ===== end NEW BLOCK ========================================================

  myOrder = order(chrTab$MaxCount,decreasing = T)
  outTab = data.frame(matrix(nrow = 0,ncol=length(chrTab)))
  names(outTab)=names(chrTab)
  outBED = data.frame(matrix(nrow = 0,ncol=length(chrBed)))
  names(outBED)=names(chrBed)
  outBinned = list()
  seen = 0
  for (i in myOrder)
  {
    #print(paste("i",i,"nrow",nrow(outTab)))
    seen = seen + 1
    if (seen %% 1000 ==0)
    {
      print(paste(seen,"/",nrow(chrTab),legit,tooClose,length(outBinned)))
    }
    if (is.na(chrTab[i,"TopBarcode"]))
      next

    # ===== NEW: drop positions not on a TTAA. Placed BEFORE the binning block =
    # so a non-TTAA position is neither kept NOR added as an anchor (outBinned),
    # i.e. it cannot suppress a nearby genuine TTAA position. Skipping i also
    # drops the matching chrBed[i,] row -> poses and bed stay row-aligned.
    if (ttaaFilter && !ttaaOK[i]) {
      nonTTAA = nonTTAA + 1
      next
    }
    # ===== end NEW ============================================================

    isTooClose = F
    myBin = round(chrTab[i,"V2"]/100)
    for (otherBin in c(myBin-1,myBin,myBin+1))
    {
      if (!is.null(outBinned[[as.character(otherBin)]]))
      {
        mySeen = outBinned[[as.character(otherBin)]]
        for (j in 1:nrow(mySeen))
        {
          #print(paste("j",j))
          if (j>0 && abs(chrTab[i,"V2"]-mySeen[j,"V2"])<minDist && (chrTab[i,"MaxCount"]<=(mySeen[j,"MaxCount"]*minFrac)) && StrDist(chrTab[i,"TopBarcode"],mySeen[j,"TopBarcode"])<=3)
          {
            isTooClose = T
            break
          }
        }
      }
      if (isTooClose)
        break
    }
    #if (chrTab[i,"TopBarcode"]=="CCCTGACTGCTCTCAACT")
    #{
    #  print(chrTab[i,])
    #  print(isTooClose)
    #}

    if (isTooClose)
    {
      tooClose = tooClose+1
    } else {
      outBinned[[as.character(myBin)]]=addToBin(outBinned[[as.character(myBin)]],chrTab[i,])
      outBinned[[as.character(myBin-1)]]=addToBin(outBinned[[as.character(myBin-1)]],chrTab[i,])
      outBinned[[as.character(myBin+1)]]=addToBin(outBinned[[as.character(myBin+1)]],chrTab[i,])
      outTab[nrow(outTab)+1,]=chrTab[i,]
      outBED[nrow(outBED)+1,]=chrBed[i,]
      legit = legit+1
    }
  }
  #print(paste(length(allOutTab),length(outTab[order(outTab[[2]]),])))
  allOutTab=rbind(allOutTab,outTab[order(outTab[[2]]),])
  #print(paste(length(allOutBED),length(outBED[order(outBED[[2]]),])))
  allOutBED=rbind(allOutBED,outBED[order(outBED[[2]]),])
  print(paste(chr,legit,tooClose))
}
print(paste("Total",legit,"legit",tooClose,"too close"))
if (ttaaFilter) print(paste("TTAA filter: dropped",nonTTAA,"positions as non-TTAA"))   # ===== NEW
write.table(allOutTab,paste0(path,name,".filtered_barcodes_poses.txt"),col.names=F,row.names=F,quote=F,sep="\t")
write.table(allOutBED,paste0(path,name,".filtered_barcodes.bed"),col.names=F,row.names=F,quote=F,sep="\t")
