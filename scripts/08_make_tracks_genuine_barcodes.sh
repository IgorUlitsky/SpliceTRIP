#!/bin/bash
# =============================================================================
# Stage 08 — intersect each construct's genuine-barcode positions against reference
# annotation tracks (MANE genes, MANE introns, gencode noncoding, repeats)
# =============================================================================
# Produces the per-construct "<name>.<track>.txt" files that stage 10
# (table_annotation) reads. One bedtools intersect per (construct x reference track x
# strand-orientation) combination, for every construct name listed in the name-list file.
#
# Usage
#   bash 08_make_tracks_genuine_barcodes.sh <name_list.txt> <mane_bed> <mane_introns_bed> <gencode_noncoding_bed> <repeats_bed>
#
# name_list.txt: one construct name per line, e.g. "5ss.txt" (matches stage 06's <mapName>.txt
# naming -- each line's "<name>.retained_barcodes_genuine.bed" must already exist, stage 06
# output). This is igoru's make_tracks_genuine_barcodes.sh, with the 4 previously-hardcoded
# annotation bed paths (all pointing at igoru's/lidiag's home dirs) turned into required
# arguments so the config drives them instead.
# =============================================================================
set -euo pipefail

FILELIST=${1:?"name list file required"}
MANE_BED=${2:?"MANE bed required"}
MANE_INTRONS_BED=${3:?"MANE introns bed required"}
GENCODE_NONCODING_BED=${4:?"gencode noncoding bed required"}
REPEATS_BED=${5:?"repeats bed required"}

echo "File with exp.names: {$FILELIST}"

echo "Getting MANE_all tracks"
while read -r p; do bedtools intersect -s -wo -a "$p.retained_barcodes_genuine.bed" -b "$MANE_BED" > "$p.MANE_all.txt"; done < "$FILELIST"

echo "Getting MANE_exonic tracks"
while read -r p; do bedtools intersect -s -wo -a "$p.retained_barcodes_genuine.bed" -b "$MANE_BED" -split > "$p.MANE_exonic.txt"; done < "$FILELIST"

echo "Getting MANE_exonic_antisense tracks"
while read -r p; do bedtools intersect -S -wo -a "$p.retained_barcodes_genuine.bed" -b "$MANE_BED" -split > "$p.MANE_exonic_antisense.txt"; done < "$FILELIST"

echo "Getting MANE_all_antisense tracks"
while read -r p; do bedtools intersect -S -wo -a "$p.retained_barcodes_genuine.bed" -b "$MANE_BED" > "$p.MANE_all_antisense.txt"; done < "$FILELIST"

echo "Getting MANE_introns_antisense tracks"
while read -r p; do bedtools intersect -S -wo -a "$p.retained_barcodes_genuine.bed" -b "$MANE_INTRONS_BED" > "$p.MANE_introns_antisense.txt"; done < "$FILELIST"

echo "Getting MANE_introns_sense tracks"
while read -r p; do bedtools intersect -s -wo -a "$p.retained_barcodes_genuine.bed" -b "$MANE_INTRONS_BED" > "$p.MANE_introns_sense.txt"; done < "$FILELIST"

echo "Getting gencode_noncoding_all tracks"
while read -r p; do bedtools intersect -s -wo -a "$p.retained_barcodes_genuine.bed" -b "$GENCODE_NONCODING_BED" > "$p.gencode_noncoding_all.txt"; done < "$FILELIST"

echo "Getting repeats tracks"
while read -r p; do bedtools intersect -s -wo -a "$p.retained_barcodes_genuine.bed" -b "$REPEATS_BED" > "$p.repeats.txt"; done < "$FILELIST"
