#!/bin/bash
# =============================================================================
# Stage 09 — GC content of each genuine barcode's +/-500bp flanking region
# =============================================================================
# Usage
#   bash 09_gc_content_single_genuine_barcodes.sh <name> <chrom_sizes> <genome_fasta>
#
# Inputs   <name>.retained_barcodes_genuine.bed   (stage 06 output)
# Outputs  <name>_with_gc.bed   (stage 06 bed + GC%% column; read by stage 10)
#
# Run once PER construct (<name> = "<mapName>.txt", matching stage 06's naming) --
# unlike stage 08/10 which take a whole name-list file in one call, this one processes a
# single construct per invocation, matching the original per-construct bsub usage.
#
# This is igoru's gc_content_single_genuine_barcodes.sh, with the 2 previously-hardcoded
# reference paths (chrom sizes, genome fasta) turned into required arguments so the config
# drives them instead. `ml BEDTools` intentionally dropped here (unlike the original) --
# this repo's convention is that required modules are loaded once up front by the caller
# (see README "Running"), not re-loaded per stage script.
# =============================================================================
set -euo pipefail

NAME=${1:?"construct name required"}
CHROM_SIZES=${2:?"chrom sizes file required"}
GENOME_FASTA=${3:?"genome fasta required"}

echo "Processing file: $NAME"
echo "Increasing Coordinates by 500 bp Upstream and Downstream"
bedtools slop -b 500 -g "$CHROM_SIZES" -i "$NAME.retained_barcodes_genuine.bed" > "${NAME}_tmp_extended.bed"

echo "Extracting FASTA"
bedtools getfasta -fi "$GENOME_FASTA" -bed "${NAME}_tmp_extended.bed" -name | \
  awk '/^>/ {sub(/::.*/, "", $0); print; next} {print}' > "${NAME}_sequences.fasta"

rm -f "${NAME}_gc_content.txt"

echo "Calculating GC Content"
awk '
/^>/ { name=substr($0,2); next }
{ seq=$0; total=length(seq); gc=gsub(/[GCgc]/,"",seq)
  printf "%s\t%.2f\n", name, (total>0 ? gc/total*100 : 0) }
' "${NAME}_sequences.fasta" > "${NAME}_gc_content.txt"

paste "$NAME.retained_barcodes_genuine.bed" "${NAME}_gc_content.txt" > "${NAME}_with_gc.bed"

echo "Checking discrepancies between the 4th and 7th columns"
awk '$4 != $7' "${NAME}_with_gc.bed"

echo "Done processing $NAME"
echo "--------------------------"
