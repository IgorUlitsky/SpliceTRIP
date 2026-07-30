library(DescTools)

path="/home/labs/ulitsky/shared/data/human/TRIP/IntronTRIP_combined/"
nameMap = "A1.txt"
barcodeLen = 18

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  nameMap = myArgs[2]
  # myArgs[3] is optional (barcodeLen keeps its default of 18 above when omitted) -- indexing
  # a shorter myArgs vector out of bounds returns NA rather than erroring, and
  # as.integer(NA) is silently NA, not 18. That NA then propagates through
  # nchar(barcodes)==barcodeLen in filterBarcodes() as an all-NA comparison; barcodes[NA]
  # keeps its original length (just filled with NA) so the function's own length==0 guard
  # never trips, but the downstream tapply() drops the NA-keyed group entirely, silently
  # producing zero rows for every position and an empty final output with no error at all
  # (confirmed: this is exactly what happened running with only 2 args, path and nameMap).
  if (length(myArgs)>=3) barcodeLen = as.integer(myArgs[3])
}
print(path)
print(nameMap)

# Read and merge all replicate bed files listed in nameMap
print('Reading data')
datas = list()
for (fileName in read.table(paste0(path, nameMap))$V1)
{
  if (file.size(paste0(path,fileName))>0)
  {
    data <- read.table(paste0(path,fileName), sep = "\t", stringsAsFactors = FALSE)
    colnames(data) <- c('chr', 'start', 'end', 'barcode', 'count', 'strand')
    datas[[fileName]]=data
  }
}
bound = do.call(rbind,datas)

# Pool all barcode:count calls at each genomic position across replicates
aggregated = aggregate(bound$barcode,by=bound[,c(1:3,6)],FUN=paste0,collapse=",")
aggregated$count=sapply(aggregated$x,function(a){length(strsplit(a,",")[[1]])})
aggregated=aggregated[,c("chr","start","end","x","count","strand")]
names(aggregated) <- c('chr', 'start', 'end', 'barcode', 'count', 'strand')

# For each position: keep only well-formed barcodes (length barcodeLen, no N),
# aggregate counts per unique barcode, and merge near-duplicate barcodes (edit distance
# <= 2) into the higher-count variant.
# Uses plain vectors + tapply instead of aggregate()/data.frame() per position
# -- with ~140k positions, aggregate()'s generic dispatch overhead dominated
# runtime (see profiling); this keeps the same grouping/merge logic (including
# its existing quirk of merging into a barcode that later gets dropped itself)
# while avoiding that per-call overhead.
filterBarcodes <- function(chr,start,end,strand,myStr)
{
  pairs = strsplit(strsplit(myStr,",")[[1]],":")
  pairs = pairs[lengths(pairs)==2]
  if (length(pairs)==0) return(NULL)
  barcodes = vapply(pairs,`[`,character(1),1)
  counts = as.numeric(vapply(pairs,`[`,character(1),2))
  ok = nchar(barcodes)==barcodeLen & !grepl("N",barcodes)
  barcodes = barcodes[ok]
  counts = counts[ok]
  if (length(barcodes)==0) return(NULL)

  sums = tapply(counts,barcodes,sum)
  ord = order(sums,decreasing=T)
  ubc = names(sums)[ord]
  usum = as.numeric(sums)[ord]
  if (length(ubc)>1)
  {
    keep = rep(T,length(ubc))
    for (i in 2:length(ubc))
    {
      for (j in 1:(i-1))
      {
        if (StrDist(ubc[i],ubc[j])<=2)
        {
          usum[j]=usum[i]+usum[j]
          keep[i]=F
          break;
        }
      }
    }
    ubc=ubc[keep]
    usum=usum[keep]
  }
  list(chr=rep(chr,length(ubc)),start=rep(start,length(ubc)),end=rep(end,length(ubc)),
       Group.1=ubc,x=usum,strand=rep(strand,length(ubc)))
}

filteredBarcodes = lapply(seq_len(nrow(aggregated)), function(i)
  filterBarcodes(aggregated$chr[i],aggregated$start[i],aggregated$end[i],aggregated$strand[i],aggregated$barcode[i]))
filteredBarcodes = filteredBarcodes[!vapply(filteredBarcodes,is.null,logical(1))]
filtered = data.frame(
  chr = unlist(lapply(filteredBarcodes,`[[`,"chr")),
  start = unlist(lapply(filteredBarcodes,`[[`,"start")),
  end = unlist(lapply(filteredBarcodes,`[[`,"end")),
  Group.1 = unlist(lapply(filteredBarcodes,`[[`,"Group.1")),
  x = unlist(lapply(filteredBarcodes,`[[`,"x")),
  strand = unlist(lapply(filteredBarcodes,`[[`,"strand")),
  stringsAsFactors = FALSE
)

# Split into barcodes with a single genomic integration site (singletons)
# vs. barcodes that map to multiple sites (multi-mappers)
locationCounts = table(filtered$Group.1)

singletons = subset(filtered,locationCounts[filtered$Group.1]==1)
write.table(singletons[order(singletons$chr,as.numeric(singletons$start)),], paste0(path,nameMap,".retained_barcodes_genuine.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
forUCSC = singletons[order(singletons$chr,as.numeric(singletons$start)),]
forUCSC[,5]=sapply(forUCSC[,5],min,1000)
write.table(forUCSC, paste0(path,nameMap,".retained_barcodes_genuine_forUCSC.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

multiMappers = subset(filtered,locationCounts[filtered$Group.1]>1)
write.table(multiMappers[order(multiMappers$chr,as.numeric(multiMappers$start)),], paste0(path,nameMap,".retained_barcodes_multimappers.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
forUCSC = multiMappers[order(multiMappers$chr,as.numeric(multiMappers$start)),]
forUCSC[,5]=sapply(forUCSC[,5],min,1000)
write.table(forUCSC, paste0(path,nameMap,".retained_barcodes_multimappers_forUCSC.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

print('Saving results')
