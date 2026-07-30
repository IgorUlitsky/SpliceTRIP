library(DescTools)
# NOTE: Biostrings is only needed for reverse-complement; it is loaded lazily inside the
# doReverseComplement block below, so the mismatch path runs with base R only.

names_with_P_list = c('5ss_Expression_rep1_S56', '5ss_Expression_rep2_S57',
                      'CPA_Expression_rep2_S54', 'CPA_Expression_rep3_S55')

buildBED <- function(df) {
  bedDF = data.frame(chr = df$Chr, start = df$Pos - 1, end = df$Pos + 1, name = rownames(df),
                       score = df$exp,
                       strand = df$Strand)
  bedDF[order(df$Chr, df$Pos), ]
}

effLength=16
myArgs = commandArgs(trailingOnly = T)
ind=1
if (length(myArgs)>0)
{
  path = myArgs[1]
  path_norm = myArgs[2]
  nameMap = myArgs[3]
  nameNorm = myArgs[4]
  nameExp = myArgs[5]
  effLength = myArgs[6]
  doReverseComplement = myArgs[7]
}
#path="/home/labs/ulitsky/shared/data/human/TRIP/set1_combined_1//"
#path_norm="/home/labs/ulitsky/shared/data/human/TRIP/set1_combined_1//"
#nameMap="CPA_noP.txt"
#nameNorm="Norm_5ss_noP_pool1_rep2_S3"
#nameExp="3intron_CPA_noP_Pool1_Rep2_S104"
#effLength=18
#doReverseComplement=F

# ============================================================================
# ===== NEW BLOCK (MM: 1-nt mismatch-tolerant matching) — optional flag =======
# Usage:
#   original behaviour (EXACT match) — 7 args:
#     Rscript normalizeTRIP_v2_MM.R <path> <path_norm> <nameMap> <nameNorm> <nameExp> <effLength> <doReverseComplement>
#   allow ONE-nt mismatch — add the literal flag "oneMismatch" as the 8th argument:
#     Rscript normalizeTRIP_v2_MM.R <...same 7 args...> oneMismatch
#
# The tolerance is HARD-FIXED at exactly ONE substitution (Hamming distance 1). This flag
# only switches that on/off; there is deliberately NO way to request 2+ mismatches. Without
# the flag the script behaves IDENTICALLY to the original normalizeTRIP_v2.R.
oneMismatch = FALSE
if (length(myArgs) >= 8 && myArgs[8] == "oneMismatch") oneMismatch = TRUE
print(paste("oneMismatch flag:", oneMismatch))
# ===== end NEW BLOCK ========================================================
# ============================================================================

# ============================================================================
# ===== NEW BLOCK (MM) — neighbor-hash helpers (FIXED at 1 mismatch) =========
# These functions are self-contained (defined here, nothing external is sourced). The
# neighbor-hash idea is adapted from
#   igoru/R/trip/cDNA/analyzeCDNAStats_run2.R  (make_neighbor_mapping / find_closest_kosher).
# Like that source, the distance is fixed at ONE substitution. Differences from it:
#   (a) ambiguous ties (a read 1 nt from >1 genuine) are DROPPED, not "first wins";
#   (b) exact matches (distance 0) always take priority over 1-nt ones;
#   (c) guarded against empty input.
# Distance = Hamming (substitutions only); all barcodes here are exactly effLength long.

# The barcode itself (distance 0) + all single-substitution variants (distance 1).
mmVariants1 <- function(bc, alphabet = c("A","C","G","T")) {
  ch = strsplit(bc, "", fixed = TRUE)[[1]]; L = length(ch)
  vars = bc; dists = 0L                                   # distance 0 = the barcode itself
  for (p in 1:L) {
    for (a in alphabet[alphabet != ch[p]]) {              # 3 alternative bases at position p
      nc = ch; nc[p] = a
      vars = c(vars, paste(nc, collapse = "")); dists = c(dists, 1L)
    }
  }
  data.frame(variant = vars, dist = dists, genuine = bc, stringsAsFactors = FALSE)
}

# Build the variant -> genuine lookup. Exact (distance 0) wins; a variant that is 1 nt
# from MORE THAN ONE genuine is ambiguous and excluded. Named vector (names = variant).
mmBuildVariantMap <- function(genuine) {
  if (length(genuine) == 0) return(character(0))
  allv = do.call(rbind, lapply(genuine, mmVariants1))
  allv = allv[order(allv$variant, allv$dist), ]          # within each variant, distance ascending
  grp = split(seq_len(nrow(allv)), allv$variant)
  res = vapply(grp, function(ix) {
    dmin = allv$dist[ix][1]                               # first = minimum distance (0 or 1)
    g = unique(allv$genuine[ix][allv$dist[ix] == dmin])  # genuine(s) at that minimum distance
    if (length(g) == 1) g else NA_character_             # exactly one -> assign; tie -> drop
  }, character(1))
  res[!is.na(res)]
}

# Look up barcodes in the variant map (NA if no match within 1 nt, or ambiguous).
mmSnap <- function(barcodes, variantMap) variantMap[barcodes]

# Re-key a count table (rownames = barcode, count column "x") onto genuine barcodes,
# summing the counts of all reads that snap to the same genuine barcode.
mmRemap <- function(tab, variantMap) {
  snapped = mmSnap(rownames(tab), variantMap)
  keep = !is.na(snapped)
  tab = tab[keep, , drop = FALSE]; snapped = snapped[keep]
  if (nrow(tab) == 0) return(tab)
  agg = aggregate(tab$x, list(snapped), FUN = sum)
  rownames(agg) = agg$Group.1; agg$Group.1 = NULL
  agg
}
# ===== end NEW BLOCK ========================================================
# ============================================================================

print(nameNorm)
print(nameExp)
print(nameMap)

mapTab=read.table(paste0(path,nameMap,".retained_barcodes_genuine.bed"),sep="\t")
name = nameExp

mapTab$V4 = substr(mapTab$V4,1,effLength)
mapTab$Ind = 1:nrow(mapTab)
barcodeToInd = aggregate(mapTab$Ind,by=list(mapTab$V4),FUN=c)
rownames(barcodeToInd)=barcodeToInd$Group.1
head(barcodeToInd)

# ===== NEW BLOCK (MM) — build the variant->genuine map from the genuine map barcodes =====
variantMap = NULL
if (oneMismatch) {
  variantMap = mmBuildVariantMap(rownames(barcodeToInd))
  print(paste("MM: built variant map with", length(variantMap),
              "unambiguous keys from", nrow(barcodeToInd), "genuine barcodes"))
}
# ===== end NEW BLOCK ====================================================================

# Filtering minNormReads = 2
minNormReads = 2

normTab=read.table(paste0(path_norm,nameNorm,".counts.txt.gz"),sep="\t",row.names = 1)
print(paste(format(nrow(normTab),big.mark = ","),"normalization barcodes"))
table(nchar(rownames(normTab)))
normTab=subset(normTab,nchar(rownames(normTab))>=effLength & nchar(rownames(normTab))<=effLength)
normTab=subset(normTab,normTab$V2>=minNormReads)
print(paste(format(nrow(normTab),big.mark=",")," normalized barcodes pre filtering"))

normTab = aggregate(normTab$V2,list(substr(rownames(normTab),1,effLength)),FUN=sum)
rownames(normTab)=normTab$Group.1
normTab$Group.1=NULL

# ===== NEW BLOCK (MM) — snap Norm barcodes onto genuine (rescue 1-nt typos, sum counts) =====
if (oneMismatch) {
  before = nrow(normTab)
  normTab = mmRemap(normTab, variantMap)
  print(paste("MM: Norm barcodes", before, "->", nrow(normTab), "after snapping to genuine"))
}
# ===== end NEW BLOCK ======================================================================

normTab=subset(normTab,rownames(normTab) %in% rownames(barcodeToInd))

print(paste(format(nrow(normTab),big.mark=",")," normalized barcodes after retaining mapped ones"))

print("Loading expression")

# Aggregating the barcodes by the first effLegnth bases
expTab=read.table(paste0(path_norm,nameExp,".counts.txt.gz"),sep="\t",row.names = 1)
expTab = aggregate(expTab$V2,list(substr(rownames(expTab),1,effLength)),FUN=sum)
rownames(expTab)=expTab$Group.1
expTab$Group.1=NULL

expTab=subset(expTab,nchar(rownames(expTab))==effLength)
print(paste(format(nrow(expTab),big.mark=",")," expression barcodes"))

if (doReverseComplement)
{
  library(Biostrings)                                    # ===== NEW: loaded lazily (only if needed)
  rownames(expTab)=as.character(reverseComplement(DNAStringSet(rownames(expTab))))
}

# ===== NEW BLOCK (MM) — snap Exp barcodes onto genuine (after any reverse-complement) =====
if (oneMismatch) {
  before = nrow(expTab)
  expTab = mmRemap(expTab, variantMap)
  print(paste("MM: Exp barcodes", before, "->", nrow(expTab), "after snapping to genuine"))
}
# ===== end NEW BLOCK =====================================================================

paste(format(length(intersect(rownames(expTab),rownames(normTab))),big.mark = ","),"barcodes in common between expression and normalization")
paste(format(length(intersect(rownames(expTab),rownames(barcodeToInd))),big.mark = ","),"barcodes in common between expression and mapping")
paste(format(length(intersect(rownames(normTab),rownames(barcodeToInd))),big.mark = ","),"barcodes in common between normalization and mapping")


print("Table of barcodeToInd indices")
print(table(sapply(barcodeToInd[rownames(normTab),"x"],length)))

# Retain single-position barcodes
sharedBarcodes=rownames(normTab)[sapply(barcodeToInd[rownames(normTab),"x"],length)==1]
length(sharedBarcodes)
df = data.frame(row.names=sharedBarcodes)
df$norm = normTab[sharedBarcodes,"x"]
df$exp = expTab[sharedBarcodes,"x"]
mapsExtracted = mapTab[unlist(barcodeToInd[match(sharedBarcodes,rownames(barcodeToInd)),"x"]),]
dim(df)
dim(mapsExtracted)
df$Chr  = mapsExtracted$V1
df$Pos  = mapsExtracted$V2
df$Strand  = mapsExtracted$V6
df$map  = mapsExtracted$V5
df[is.na(df)]=0
print(paste("Writing table with ",nrow(df)," entries"))
write.table(buildBED(df), paste0(path, name, ".bed"), quote = F, col.names = F, row.names = F, sep = "\t")
write.table(df,paste0(path,name,".combined_counts.txt"),sep="\t",quote=F)
