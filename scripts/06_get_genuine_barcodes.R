library(DescTools)

path="/home/labs/ulitsky/shared/data/human/TRIP/IntronTRIP_combined/"
nameMap = "A1.txt"
barcodeLen = 18
dpniiBedPath = "/home/projects/ulitsky/igoru/data/human/hg38_dpnii.bed"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  nameMap = myArgs[2]
  # myArgs[3]/[4] are optional (barcodeLen/dpniiBedPath keep their defaults above when
  # omitted) -- indexing a shorter myArgs vector out of bounds returns NA rather than
  # erroring, and as.integer(NA) is silently NA, not 18. That NA then propagates through
  # nchar(barcodes)==barcodeLen in filterBarcodes() as an all-NA comparison; barcodes[NA]
  # keeps its original length (just filled with NA) so the function's own length==0 guard
  # never trips, but the downstream tapply() drops the NA-keyed group entirely, silently
  # producing zero rows for every position and an empty final output with no error at all
  # (confirmed: this is exactly what happened running with only 2 args, path and nameMap).
  if (length(myArgs)>=3) barcodeLen = as.integer(myArgs[3])
  if (length(myArgs)>=4 && nzchar(myArgs[4])) dpniiBedPath = myArgs[4]
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

# ============================================================================
# ===== NEW BLOCK (DpnII distance) — for each genuine (singleton) position, distance
# to the nearest DOWNSTREAM DpnII site, downstream relative to the integration's own
# strand (+ strand: next site at a higher coordinate; - strand: next site at a lower
# coordinate). hg38_dpnii.bed lists DpnII RESTRICTION FRAGMENTS (a Hi-C-style digest
# map: "chr1 0 11159 HIC_chr1_1 ..."), not individual cut-site coordinates directly --
# but since consecutive fragments share a boundary and that shared boundary IS the cut
# site, the site positions are just every fragment's start EXCEPT each chromosome's
# very first fragment (whose start is 0, a chromosome-boundary artifact, not a real
# DpnII site).
print(paste("Loading DpnII sites from", dpniiBedPath))
dpnii = read.table(dpniiBedPath, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = c("character","integer","NULL","NULL","NULL","NULL"))
colnames(dpnii) = c("chrom","start")
dpnii = dpnii[order(dpnii$chrom, dpnii$start), ]
dpnii = dpnii[duplicated(dpnii$chrom), ]                # drop each chrom's first (boundary) fragment
dpniiByChrom = lapply(split(dpnii$start, dpnii$chrom), sort)
print(paste("Loaded", nrow(dpnii), "DpnII sites on", length(dpniiByChrom), "contigs"))

# vec is sorted ascending 0-based DpnII site starts for one contig (NULL if the contig
# isn't in the reference at all). NA if pos is past the last / before the first site.
nearestDownstreamDpnII <- function(chrom, pos, strand, byChrom)
{
  vec = byChrom[[chrom]]
  if (is.null(vec)) return(NA_real_)
  if (strand == "+") {
    idx = findInterval(pos, vec) + 1L
    if (idx > length(vec)) return(NA_real_)
    vec[idx] - pos
  } else {
    idx = findInterval(pos - 1L, vec)
    if (idx < 1L) return(NA_real_)
    pos - vec[idx]
  }
}

singletonsOrdered = singletons[order(singletons$chr,as.numeric(singletons$start)),]
dpniiDist = mapply(nearestDownstreamDpnII, singletonsOrdered$chr, as.numeric(singletonsOrdered$start),
                    singletonsOrdered$strand, MoreArgs = list(byChrom = dpniiByChrom))

# A handful of distances come out implausibly large (multi-kb to tens of kb) -- almost
# certainly assembly gaps (N stretches) or other low-complexity/repeat-masked regions
# where the reference genuinely has no GATC for a long stretch, not a real biological
# signal. Cap at 10x the by-chance expected distance (assuming a gap-free genome at
# genomeGC% GC content, independent bases: P(GATC) = (GC/2)^2 * ((1-GC)/2)^2, expected
# distance = 1/P(GATC)) and set anything beyond that to NA, excluded from downstream
# means/plots rather than distorting them with a few genome-gap outliers.
genomeGC = 0.44
pGATC = (genomeGC/2)^2 * ((1-genomeGC)/2)^2
expectedDist = 1/pGATC
maxPlausibleDist = 10 * expectedDist
nCapped = sum(!is.na(dpniiDist) & dpniiDist > maxPlausibleDist)
dpniiDist[!is.na(dpniiDist) & dpniiDist > maxPlausibleDist] = NA_real_
print(paste0("Expected by-chance DpnII spacing at ", genomeGC*100, "% GC (no Ns): ",
             round(expectedDist,1), " bp; NA-ing ", nCapped,
             " distances beyond 10x that (", round(maxPlausibleDist,1), " bp)"))

print(paste("DpnII downstream distance: median", median(dpniiDist, na.rm = TRUE),
            "bp across", sum(!is.na(dpniiDist)), "of", length(dpniiDist), "positions"))
write.table(
  data.frame(chr = singletonsOrdered$chr, start = singletonsOrdered$start,
             strand = singletonsOrdered$strand, barcode = singletonsOrdered$Group.1,
             dpnII_dist = dpniiDist),
  paste0(path, nameMap, ".retained_barcodes_genuine.dnpII_dists.txt"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# ===== end NEW BLOCK ========================================================
# ============================================================================

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
