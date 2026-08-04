#####
# Regenerates the MANE introns bed file directly from MANE.GRCh38.v1.0.refseq_genomic.bed's
# own BED12 exon block structure, so it fits that file exactly (same transcript set,
# correct 0-based coordinates) -- replacing
# /home/projects/ulitsky/lidiag/genomic_data/GCF_000001405.40_GRCh38.p14_MANE_introns_new_format.bed,
# which has a baked-in V2 (start) 1-based/0-based mix-up (see
# 11_table_annotation_genuine_barcodes.R's earlier fix) and covers a stale transcript set:
# 136 transcripts not in the current maneBED at all, and 79 multi-exon maneBED transcripts
# missing entirely (versus 1025 single-exon maneBED transcripts correctly absent, since
# those have no introns).
#
# Format matches the original file exactly: chrom, start, end, "<transcript>.<intronNum>",
# score(=intronNum), strand -- intron numbering confirmed empirically to run in genomic
# (increasing-coordinate) order regardless of strand, matching BED12's own block order.
#
# Usage: Rscript make_mane_introns.R [mane_bed] [out_bed]
# Both args optional, defaulting to the paths below.
#####
library(data.table)

maneBedPath = "/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed"
outFile = "/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.introns_new_format.bed"

myArgs = commandArgs(trailingOnly = TRUE)
if (length(myArgs) >= 1 && nzchar(myArgs[1])) maneBedPath = myArgs[1]
if (length(myArgs) >= 2 && nzchar(myArgs[2])) outFile = myArgs[2]

maneBED = fread(maneBedPath, sep = "\t", header = FALSE)

multiExon = maneBED[V10 > 1]
cat(sprintf("%d/%d MANE transcripts have >1 exon (will get introns)\n", nrow(multiExon), nrow(maneBED)))

deriveIntrons = function(chrom, chromStart, strand, txName, blockSizesStr, blockStartsStr)
{
  sizes = as.integer(strsplit(blockSizesStr, ",", fixed = TRUE)[[1]])
  starts = as.integer(strsplit(blockStartsStr, ",", fixed = TRUE)[[1]])
  exonStart = chromStart + starts
  exonEnd = exonStart + sizes
  n = length(exonStart)
  intronStart = exonEnd[1:(n - 1)]
  intronEnd = exonStart[2:n]
  data.table(chrom = chrom, start = intronStart, end = intronEnd,
             name = paste0(txName, ".", seq_len(n - 1)),
             score = seq_len(n - 1), strand = strand)
}

introns = rbindlist(lapply(seq_len(nrow(multiExon)), function(i) {
  r = multiExon[i]
  deriveIntrons(r$V1, r$V2, r$V6, r$V4, r$V11, r$V12)
}))

setorder(introns, chrom, start)

cat(sprintf("%d introns derived across %d transcripts\n", nrow(introns), length(unique(multiExon$V4))))

write.table(introns, outFile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
cat(sprintf("wrote %s\n", outFile))
