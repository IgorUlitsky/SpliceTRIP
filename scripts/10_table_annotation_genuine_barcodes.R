library(rtracklayer)
library(dplyr)

# Path to dir with TRIP data, and to the name-list file (one construct name per line, e.g.
# "5ss.txt" -- matches stage 06's <mapName>.txt naming and stage 08's name-list file).
path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined//"
namesMaps_file <- "names.complete_maps.txt"

# Reference annotation paths -- all previously hardcoded to igoru's/lidiag's home dirs;
# now optional CLI args (in this fixed order) so a config file drives them, defaulting to
# the original paths when omitted so this still runs unmodified from the command line.
maneBedPath = "/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed"
maneIntronsBedPath = "/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.introns_new_format.bed"
gencodeNoncodingBedPath = "/home/projects/ulitsky/lidiag/genomic_data/gencode.v46.annotation.noncoding.bed"
trans2genePath = "/home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.trans2gene.txt"
maneExpressionPath = "/home/projects/ulitsky/igoru/data/human/ENCODE/K562/K562_1_MANEReadsPerGene.out.tab"
chromHmmBedPath = "/home/projects/ulitsky/shared/data/human/K562_ChromHMM_Kellis_ENCFF319VXX.bed.gz"
spinBedPath = "/home/projects/ulitsky/igoru/R/trip/data/K562_SPIN_25kb_hg38.bed"

# Arguments from the terminal
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  namesMaps_file = myArgs[2]
  # Args 3-9 are all optional (each keeps its default above when omitted) -- same
  # length-guarded pattern as stage 06's barcodeLen/dpniiBedPath, so a missing trailing
  # arg falls back to the documented default instead of silently becoming NA (as.integer's
  # or as.character's NA on an out-of-bounds myArgs[i] is silent, not an error).
  if (length(myArgs)>=3 && nzchar(myArgs[3])) maneBedPath = myArgs[3]
  if (length(myArgs)>=4 && nzchar(myArgs[4])) maneIntronsBedPath = myArgs[4]
  if (length(myArgs)>=5 && nzchar(myArgs[5])) gencodeNoncodingBedPath = myArgs[5]
  if (length(myArgs)>=6 && nzchar(myArgs[6])) trans2genePath = myArgs[6]
  if (length(myArgs)>=7 && nzchar(myArgs[7])) maneExpressionPath = myArgs[7]
  if (length(myArgs)>=8 && nzchar(myArgs[8])) chromHmmBedPath = myArgs[8]
  if (length(myArgs)>=9 && nzchar(myArgs[9])) spinBedPath = myArgs[9]
}
namesMaps<-read.table(paste0(path, namesMaps_file))$V1
print(path)
print(namesMaps)

# External annotation data
trans2gene = read.table(trans2genePath,sep="\t")
maneExp = read.table(maneExpressionPath,sep="\t",header=F)
maneBED = read.table(maneBedPath,sep="\t",header=F)
maneBED$Length = maneBED$V3-maneBED$V2
maneBED$ExonicLength = sapply(strsplit(as.character(maneBED$V11),","),function(x){sum(as.numeric(x))})
# MANE.GRCh38.v1.0.refseq_genomic.introns_new_format.bed is generated directly from
# maneBedPath's own BED12 exon block structure (see make_mane_introns.R) -- correct
# 0-based coordinates by construction and covers exactly maneBedPath's transcript set, so
# (unlike the old lidiag reference file this replaced) no coordinate correction is needed
# here.
intronBED = read.table(maneIntronsBedPath,sep="\t",header=F)
intronBED$Length = intronBED$V3-intronBED$V2
gencodeBED = read.table(gencodeNoncodingBedPath,sep="\t",header=F)
gencodeBED$Length = gencodeBED$V3-gencodeBED$V2

bed_chrom_hmm <- gzfile(chromHmmBedPath)
chrom_hmm = import(bed_chrom_hmm)

spin = import(spinBedPath)

AnnotTabs = list()
repeats = list()
GcTabs = list()
for (name in namesMaps)
{
  print(name)
  # Reading data
  AnnotTabs[[name]] = read.table(paste0(path, name, ".retained_barcodes_genuine.bed"), sep = "\t", header = T)
  names(AnnotTabs[[name]])=c("Chr","Pos","End","Barcode","Score","Strand")
  rownames(AnnotTabs[[name]])=AnnotTabs[[name]]$Barcode

  maneSense = read.table(paste0(path, name, ".MANE_all.txt"), sep = "\t", header = F,quote="",comment="")
  maneSense[941:943,]
  maneExonic = read.table(paste0(path, name, ".MANE_exonic.txt"), sep = "\t", header = F)
  maneAntisense = read.table(paste0(path, name, ".MANE_all_antisense.txt"), sep = "\t", header = F)
  maneExonicAntisense = read.table(paste0(path, name, ".MANE_exonic_antisense.txt"), sep = "\t", header = F)
  intronSense = read.table(paste0(path, name, ".MANE_introns_sense.txt"), sep = "\t", header = F)
  intronAntisense = read.table(paste0(path, name, ".MANE_introns_antisense.txt"), sep = "\t", header = F)
  gencode_noncoding = read.table(paste0(path, name, ".gencode_noncoding_all.txt"), sep = "\t", header = F)

  # Add MANE data to AnnotTabs
  AnnotTabs[[name]]$MANE_all = maneSense[match(rownames(AnnotTabs[[name]]), maneSense$V4), "V10"]
  AnnotTabs[[name]]$MANE_exonic = maneExonic[match(rownames(AnnotTabs[[name]]), maneExonic$V4), "V10"]
  AnnotTabs[[name]]$MANE_antisense = maneAntisense[match(rownames(AnnotTabs[[name]]), maneAntisense$V4), "V10"]
  AnnotTabs[[name]]$MANE_exonic_antisense = maneExonicAntisense[match(rownames(AnnotTabs[[name]]), maneExonicAntisense$V4), "V10"]
  AnnotTabs[[name]]$MANE_intron = intronSense[match(rownames(AnnotTabs[[name]]), intronSense$V4), "V10"]
  AnnotTabs[[name]]$MANE_intron_antisense = intronAntisense[match(rownames(AnnotTabs[[name]]), intronAntisense$V4), "V10"]
  AnnotTabs[[name]]$gencode_noncoding = gencode_noncoding[match(rownames(AnnotTabs[[name]]), gencode_noncoding$V4), "V10"]
  AnnotTabs[[name]]$annot_rna_type = gencode_noncoding[match(rownames(AnnotTabs[[name]]), gencode_noncoding$V4), "V13"]

  AnnotTabs[[name]] <- AnnotTabs[[name]] %>%
    mutate(
      annot_rna_type = as.character(annot_rna_type),
      annot_rna_type = case_when(
        is.na(annot_rna_type) & !is.na(MANE_all) ~ "MANE_mRNA",
        TRUE ~ annot_rna_type
      )
    )
  AnnotTabs[[name]]$gencode_noncoding_name = gencode_noncoding[match(rownames(AnnotTabs[[name]]), gencode_noncoding$V4), "V14"]

  # Add gene symbols and expression data
  AnnotTabs[[name]]$MANE_all_symbol = trans2gene[match(AnnotTabs[[name]]$MANE_all, trans2gene$V2), "V1"]
  AnnotTabs[[name]]$MANE_all_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_all, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  AnnotTabs[[name]]$MANE_exonic_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_exonic, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  AnnotTabs[[name]]$MANE_antisense_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_antisense, trans2gene$V2), "V1"], maneExp$V1), "V2"]
  AnnotTabs[[name]]$MANE_exonic_antisense_exp = maneExp[match(trans2gene[match(AnnotTabs[[name]]$MANE_exonic_antisense, trans2gene$V2), "V1"], maneExp$V1), "V2"]

  # Normalize expression data
  AnnotTabs[[name]]$MANE_all_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_all_exp / maneBED[match(AnnotTabs[[name]]$MANE_all, maneBED$V4), "ExonicLength"])
  AnnotTabs[[name]]$MANE_all_exp.norm.log10 <- log((AnnotTabs[[name]]$MANE_all_exp.norm + 1), 10)
  AnnotTabs[[name]]$MANE_exonic_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_exonic_exp / maneBED[match(AnnotTabs[[name]]$MANE_exonic, maneBED$V4), "ExonicLength"])
  AnnotTabs[[name]]$MANE_antisense_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_antisense_exp / maneBED[match(AnnotTabs[[name]]$MANE_antisense, maneBED$V4), "ExonicLength"])
  AnnotTabs[[name]]$MANE_exonic_antisense_exp.norm = 1000 * (AnnotTabs[[name]]$MANE_exonic_antisense_exp / maneBED[match(AnnotTabs[[name]]$MANE_exonic_antisense, maneBED$V4), "ExonicLength"])

  # Process intronic data
  AnnotTabs[[name]]$MANE_intronic = AnnotTabs[[name]]$MANE_all
  AnnotTabs[[name]]$MANE_intronic[!is.na(AnnotTabs[[name]]$MANE_exonic)] = NA

  # Calculate distances and relative positions
  AnnotTabs[[name]]$MANElength = maneBED[match(AnnotTabs[[name]]$MANE_all, maneBED$V4), "Length"]
  AnnotTabs[[name]]$MANE_intron_length = intronBED[match(AnnotTabs[[name]]$MANE_intron, intronBED$V4), "Length"]
  # Matched on MANE_all (not MANE_intronic): MANE_intronic is NA'd out for exonic
  # positions (see above), so matching on it here silently dropped DistToStart/DistToEnd
  # for any barcode landing in a MANE exon -- MANElength/MANE_all_symbol stayed correct
  # (they're matched on MANE_all), but the later gencode_noncoding fallback below would
  # then unconditionally fill the resulting NA from *any* coincidentally-overlapping
  # gencode feature, attributing that unrelated gene's distances/length under the real
  # gene's symbol (e.g. AHCTF1 shown with KIF28P's distances). Matching on MANE_all
  # directly gives every MANE-overlapping position (exonic or intronic) its own gene's
  # boundaries, so the fallback below now only fires for positions with no MANE overlap
  # at all.
  AnnotTabs[[name]]$DistToStart = (AnnotTabs[[name]]$Pos - maneBED[match(AnnotTabs[[name]]$MANE_all, maneBED$V4), "V2"])
  AnnotTabs[[name]]$DistToEnd = (maneBED[match(AnnotTabs[[name]]$MANE_all, maneBED$V4), "V3"] - AnnotTabs[[name]]$Pos)
  AnnotTabs[[name]]$DistToTSS = ifelse(AnnotTabs[[name]]$Strand == "+", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd)
  AnnotTabs[[name]]$DistToTES = ifelse(AnnotTabs[[name]]$Strand == "-", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd)
  AnnotTabs[[name]]$RelPos = AnnotTabs[[name]]$DistToTSS / AnnotTabs[[name]]$MANElength


  # Calculate distances to intron boundaries
  AnnotTabs[[name]]$DistToIntronStart = (AnnotTabs[[name]]$Pos - intronBED[match(AnnotTabs[[name]]$MANE_intron, intronBED$V4), "V2"])
  AnnotTabs[[name]]$DistToIntronEnd = (intronBED[match(AnnotTabs[[name]]$MANE_intron, intronBED$V4), "V3"] - AnnotTabs[[name]]$Pos)
  AnnotTabs[[name]]$DistToIntron5 = ifelse(AnnotTabs[[name]]$Strand == "+", AnnotTabs[[name]]$DistToIntronStart, AnnotTabs[[name]]$DistToIntronEnd)
  AnnotTabs[[name]]$DistToIntron3 = ifelse(AnnotTabs[[name]]$Strand == "-", AnnotTabs[[name]]$DistToIntronStart, AnnotTabs[[name]]$DistToIntronEnd)
  AnnotTabs[[name]]$RelIntronPos = AnnotTabs[[name]]$DistToIntron5 / AnnotTabs[[name]]$MANE_intron_length

  # Adding non-coding genes info
  AnnotTabs[[name]]$gencodelength = gencodeBED[match(AnnotTabs[[name]]$gencode_noncoding, gencodeBED$V4), "Length"]
  AnnotTabs[[name]]$DistToStart <- ifelse(
    is.na(AnnotTabs[[name]]$DistToStart),
    AnnotTabs[[name]]$Pos - gencodeBED[match(AnnotTabs[[name]]$gencode_noncoding, gencodeBED$V4), "V2"],
    AnnotTabs[[name]]$DistToStart
  )

  AnnotTabs[[name]]$DistToEnd <- ifelse(
    is.na(AnnotTabs[[name]]$DistToEnd),
    gencodeBED[match(AnnotTabs[[name]]$gencode_noncoding, gencodeBED$V4), "V3"] - AnnotTabs[[name]]$Pos,
    AnnotTabs[[name]]$DistToEnd
  )
  AnnotTabs[[name]]$DistToTSS <- ifelse(
    is.na(AnnotTabs[[name]]$DistToTSS),
    ifelse(AnnotTabs[[name]]$Strand == "+", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd),
    AnnotTabs[[name]]$DistToTSS
  )

  AnnotTabs[[name]]$DistToTES <- ifelse(
    is.na(AnnotTabs[[name]]$DistToTES),
    ifelse(AnnotTabs[[name]]$Strand == "-", AnnotTabs[[name]]$DistToStart, AnnotTabs[[name]]$DistToEnd),
    AnnotTabs[[name]]$DistToTES
  )

  AnnotTabs[[name]]$RelPos <- ifelse(
    is.na(AnnotTabs[[name]]$RelPos),
    AnnotTabs[[name]]$DistToTSS / AnnotTabs[[name]]$gencodelength,
    AnnotTabs[[name]]$RelPos
  )

  # Adding ChromHMM
  myPoses = GRanges(seqnames = AnnotTabs[[name]]$Chr, ranges = IRanges(start = AnnotTabs[[name]]$Pos, width = 1))
  overlapsChromHmm = findOverlaps(myPoses, chrom_hmm)
  AnnotTabs[[name]][queryHits(overlapsChromHmm), "ChromHMM"] = chrom_hmm[subjectHits(overlapsChromHmm)]$name

  # Adding SPIN state
  overlapsSpin = findOverlaps(myPoses, spin)
  AnnotTabs[[name]][queryHits(overlapsSpin), "SPIN"] = spin[subjectHits(overlapsSpin)]$name

  #Adding GC content
  GcTabs[[name]] <- read.table(paste0(path, name, "_with_gc.bed"), sep = "\t", header = FALSE,fill=T)
  AnnotTabs[[name]]$GC <- GcTabs[[name]]$V8[match(rownames(AnnotTabs[[name]]), GcTabs[[name]]$V7)]

  # Deviding on groups
  gc_mean <- mean(AnnotTabs[[name]]$GC, na.rm = TRUE)
  gc_sd <- sd(AnnotTabs[[name]]$GC, na.rm = TRUE)

  low_threshold <- gc_mean - gc_sd
  high_threshold <- gc_mean + gc_sd

  AnnotTabs[[name]]$GC_category <- cut(
    AnnotTabs[[name]]$GC,
    breaks = c(-Inf, low_threshold, high_threshold, Inf),
    labels = c("Low GC", "Average GC", "High GC"),
    include.lowest = TRUE
  )

  #Adding repeats
  repeats[[name]] <- read.table(paste0(path, name, ".repeats.txt"), sep = "\t", header = FALSE)
  AnnotTabs[[name]]$RepeatGroup <- repeats[[name]]$V13[match(rownames(AnnotTabs[[name]]), repeats[[name]]$V4)]

  # Write final table to file
  write.table(AnnotTabs[[name]], paste0(path, name, ".annotated_tables.txt"), quote = F, col.names = T, row.names = T, sep = "\t")
}
