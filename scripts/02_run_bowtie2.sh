#!/bin/bash
# =============================================================================
# Stage 02 — align Map libraries to hg38 and reduce read 2 to a BED
# =============================================================================
# RECONSTRUCTED SCRIPT. Unlike every other stage, no script for this step survives
# in the run directories: the alignment was run ad hoc. The commands below are the
# ones documented for the project (bowtie2 v2.2.9, paired-end, end-to-end, -X 2000,
# hg38), reassembled into a runnable script. Treat the flags as authoritative and
# the packaging as new.
#
# Usage
#   bash 02_run_bowtie2.sh <name> [threads] [bowtie2_index]
#
# Inputs   01_cutadapt/<name>.R1.fastq.gz , 01_cutadapt/<name>.R2.fastq.gz  (stage 01)
# Outputs  <name>.bam              sorted alignments of both mates
#          <name>.second.bam       read 2 only (flag 0x80 = second in pair)
#          02_mapping/<name>.second_read.bed.gz    BED12, one line per R2 alignment
#
# Only R2 is kept downstream: R1 is the constant reporter sequence (its barcode was
# already captured in stage 01), R2 is the genomic DNA next to the integration site.
# analyzeTRIPMapping (stage 03) rejoins the two by read name.
#
# Caveat worth knowing: R2 is not length-trimmed here, and bowtie2 aligns end-to-end.
# A long R2 that reads past the genomic part back into the insert therefore fails to
# align at all, which lowers the mapping rate for longer read lengths. If mapping
# percentages look low, check the R2 length / trimming before suspecting the genome.
# =============================================================================
set -euo pipefail

NAME=${1:?"sample name required"}
THREADS=${2:-10}
INDEX=${3:-/home/projects/ulitsky/shared/data/human/sequence/bowtie2/hg38}

mkdir -p 02_mapping

# 1. Align both mates, convert SAM -> BAM, sort by coordinate.
bowtie2 -p "$THREADS" -X 2000 -x "$INDEX" \
  -1 "01_cutadapt/${NAME}.R1.fastq.gz" \
  -2 "01_cutadapt/${NAME}.R2.fastq.gz" \
  2> "${NAME}_bowtie2.log" \
  | samtools view -b - \
  | samtools sort -o "${NAME}.bam" -
samtools index "${NAME}.bam"

# 2. Keep only read 2 (the genomic mate).
samtools view -f 0x80 -b "${NAME}.bam" > "${NAME}.second.bam"

# 3. One BED12 line per R2 alignment; this is what stage 03 reads.
bedtools bamtobed -bed12 -i "${NAME}.second.bam" | gzip > "02_mapping/${NAME}.second_read.bed.gz"

echo "stage 02 done: $NAME"
