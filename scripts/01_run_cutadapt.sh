#!/bin/bash
# =============================================================================
# Stage 01 — read trimming (cutadapt) + per-read barcode extraction / tallying
# =============================================================================
# In production this step exists only as per-sample `.csh` job scripts inside each
# run directory (e.g. shared/data/human/TRIP/run5/Map_5ss_noP_pool1_rep2_S166.csh,
# .../run5/count_exp.csh). This file is a faithful, argument-driven template of
# those commands - the cutadapt invocations and awk one-liners are copied verbatim
# from the real job scripts; only the sample name and the primer are parameterised.
#
# Usage
#   bash 01_run_cutadapt.sh map   <name> [threads]     # Map libraries   (DNA, paired, R2 used)
#   bash 01_run_cutadapt.sh count <name> [threads]     # Norm/Exp libraries (barcode tally only)
#
# Inputs   <name>_R1_001.fastq.gz , <name>_R2_001.fastq.gz   (in the current directory)
# Outputs  01_cutadapt/<name>.R1.fastq.gz , 01_cutadapt/<name>.R2.fastq.gz
#          map   -> <name>.barcodes.txt.gz   (read_name, UMI, 18-nt barcode)
#          count -> <name>.counts.txt.gz     (barcode, read_count; sorted descending)
#
# Library-dependent primers (override by exporting these before the call).
# The default R1 primer is the universal 28-nt one; run5 RNA-readout libraries use a
# longer 5' primer, e.g.
#   export R1_PRIMER=ACATGTCCGTAAATGCTTTCTGCCATTCCTTGAAGTTCAACTCGGTACCAGAGCCACC
# The barcode offset used when parsing the read header is derived from the primer
# length, so it follows automatically (28 -> 29, 58 -> 59).
# =============================================================================
set -euo pipefail

MODE=${1:?"mode required: map | count"}
NAME=${2:?"sample name required"}
THREADS=${3:-10}

R1_PRIMER=${R1_PRIMER:-AGAAGTTCAACTCGGTACCAGAGCCACC}      # 28 nt, universal R1 primer
R1_SUFFIX=${R1_SUFFIX:-GGATC}                             # DpnII junction (Map / IntronTRIP / run5 readouts)
R2_PRIMER=${R2_PRIMER:-GTACGTCACAATATGATTATCTTTCTAGGGTTAA} # Map only: ITR primer, ends in the piggyBac TTAA
BARCODE_LEN=${BARCODE_LEN:-18}

# 18 N's for the barcode part of the cutadapt pattern
NS=$(printf 'N%.0s' $(seq 1 "$BARCODE_LEN"))
# barcode starts right after the primer in the `barcode=` header token
OFFSET=$(( ${#R1_PRIMER} + 1 ))

mkdir -p 01_cutadapt

case "$MODE" in
map)
  # Pass 1: find UMI (--cut 4) + primer + barcode on R1; strip the ITR primer
  #         (which ends in TTAA) off R2, so R2 starts at the first genomic base.
  cutadapt -n 1 -j "$THREADS" --cut 4 \
    -g "${R1_PRIMER}${NS}${R1_SUFFIX}" \
    --rename='{id} umi={cut_prefix} adapter={adapter_name} barcode={match_sequence}' \
    -G "${R2_PRIMER}" \
    -o "./01_cutadapt/${NAME}.p1.R1.fastq.gz" -p "./01_cutadapt/${NAME}.p1.R2.fastq.gz" \
    "${NAME}_R1_001.fastq.gz" "${NAME}_R2_001.fastq.gz" \
    --discard-untrimmed > "${NAME}_cutadapt_pass1.log"

  # Pass 2: remove read-through into the construct from both mates.
  cutadapt -n 2 -j "$THREADS" -m 1:25 \
    -a TTAACCCTAGAAAGATAATCATATTGTGACGTAC \
    -A "GATCC${NS}GGTGGCTCTGGTACCGAGTTGAACTTCT" \
    -A GGTGGCTCTGGTACCGAGTTGAACTTCT \
    -o "./01_cutadapt/${NAME}.R1.fastq.gz" -p "./01_cutadapt/${NAME}.R2.fastq.gz" \
    "./01_cutadapt/${NAME}.p1.R1.fastq.gz" "./01_cutadapt/${NAME}.p1.R2.fastq.gz" \
    > "${NAME}_cutadapt_pass2.log"

  # Per-read barcode table: read_name, UMI, barcode. Only exactly-18-nt barcodes are kept.
  zcat "./01_cutadapt/${NAME}.R1.fastq.gz" \
    | sed 's/umi=//' | sed 's/barcode=//' \
    | awk -v off="$OFFSET" -v len="$BARCODE_LEN" \
        '{if (NR % 4 == 1) {print $1 "\t" $2 "\t" substr($4,off,len)}}' \
    | awk -v len="$BARCODE_LEN" 'length($3)==len' \
    | gzip > "${NAME}.barcodes.txt.gz"
  ;;

count)
  # Norm / Exp libraries: single pass, R2 is not used downstream. -g is UNANCHORED here
  # (no leading ^), matching every real historical cut_norm.csh/cut_exp.csh verified across
  # run2/run4/run7/SON_run2/IntronTRIP_run2 -- an anchored ^ was found to have crept in here
  # and would be stricter than what was actually run (requiring the primer at read position
  # 0 exactly, dropping reads with a few extra leading bases that the real commands kept).
  cutadapt -j "$THREADS" -n 1 --cut 4 \
    -g "${R1_PRIMER}${NS}${R1_SUFFIX}" \
    --rename='{id} umi={cut_prefix} adapter={adapter_name} barcode={match_sequence}' \
    -o "01_cutadapt/${NAME}.R1.fastq.gz" -p "01_cutadapt/${NAME}.R2.fastq.gz" \
    "${NAME}_R1_001.fastq.gz" "${NAME}_R2_001.fastq.gz" \
    --discard-untrimmed > "${NAME}_cutadapt.log"

  # Barcode -> read count, sorted by count descending.
  zcat "./01_cutadapt/${NAME}.R1.fastq.gz" \
    | sed 's/umi=//' | sed 's/barcode=//' \
    | awk -v off="$OFFSET" -v len="$BARCODE_LEN" '{if (NR % 4 == 1) {print substr($4,off,len)}}' \
    | sort | uniq -c | sort -k1,1rn \
    | awk '{print $2 "\t" $1}' \
    | gzip > "${NAME}.counts.txt.gz"
  ;;

*)
  echo "unknown mode '$MODE' (expected: map | count)" >&2; exit 1;;
esac

echo "stage 01 done: $NAME ($MODE)"
