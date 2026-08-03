#!/bin/bash
# One-off resume script: re-run single_combined from stage 04 onward (target=annotated),
# skipping stages 01-03 since their outputs already exist. Needed because stage 04's
# per-sample outputs currently on disk (single_combined/*.filtered_barcodes.bed, dated
# 2026-07-24) predate this repo's first commit (2026-07-30, which includes the
# position-proximity/edit-distance barcode dedup in 04's isTooClose logic) -- so
# everything downstream of stage 04 needs regenerating to pick up that fix.
#
# Unlike run_pipeline.sh's stages 01-04 (which all run from a single shared $WORK, with
# raw fastqs symlinked IN), the REAL historical data has stage 03's outputs
# (.barcodes_poses.txt/.barcodes.bed) sitting directly in each sample's own run3//run5/
# directory, not in single_combined/ -- confirmed: single_combined/ only has symlinks for
# .filtered_barcodes.bed onward, not for these. So stage 04 has to run once per FASTQ_MAPS
# sample, cd'd into that sample's own run directory (dirname of its path_prefix), and
# stages 05+ run from single_combined/ same as run_pipeline.sh (that's where the
# cross-sample manifests/symlinks/batching already live).
#
# Reuses run_pipeline.sh's own config-parsing + stage 05-10 command blocks verbatim (not a
# reimplementation) so behavior stays identical to the real pipeline; only stage 04's
# per-sample working directory differs from run_pipeline.sh's assumption.
set -euo pipefail

CONFIG=${1:?"usage: resume_from_stage04_single_combined.sh <config.txt> <single_combined_dir>"}
WORK=${2:?"usage: resume_from_stage04_single_combined.sh <config.txt> <single_combined_dir>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"

echo "=== Resume from stage 04 ==="
echo "config: $CONFIG"
echo "single_combined dir: $WORK"
echo

PARSED=$(mktemp -d)
trap 'rm -rf "$PARSED"' EXIT
for section in GLOBAL FASTQ_MAPS FASTQ_NORM FASTQ_SAMPLES COMBOS SPLICED_UNSPLICED; do
  awk -v s="### $section" '
    $0==s {grab=1; next}
    /^### /{grab=0}
    grab && NF && $0 !~ /^#/ {print}
  ' "$CONFIG" > "$PARSED/$section.tsv"
done

while IFS=$'\t' read -r key value; do
  [ -z "$key" ] && continue
  export "$key"="$value"
done < "$PARSED/GLOBAL.tsv"

: "${BARCODE_LEN:?BARCODE_LEN missing from ### GLOBAL}"
: "${CONT_FACTOR:=10}"
TTAA_BED=${TTAA_BED:-}
DPNII_BED=${DPNII_BED:-}
MANE_BED=${MANE_BED:-}
MANE_INTRONS_BED=${MANE_INTRONS_BED:-}
GENCODE_NONCODING_BED=${GENCODE_NONCODING_BED:-}
REPEATS_BED=${REPEATS_BED:-}
CHROM_SIZES=${CHROM_SIZES:-}
GENOME_FASTA=${GENOME_FASTA:-}
TRANS2GENE=${TRANS2GENE:-}
MANE_EXPRESSION=${MANE_EXPRESSION:-}
CHROMHMM_BED=${CHROMHMM_BED:-}
SPIN_BED=${SPIN_BED:-}

echo "TTAA_BED=${TTAA_BED:-<off>}  CONT_FACTOR=$CONT_FACTOR"
echo

# =============================================================================
# Stage 04: filterMappings (+ TTAA if configured) -- run per sample, in that sample's own
# run directory (dirname of its FASTQ_MAPS path_prefix), not in $WORK.
# =============================================================================
echo "=== Stage 04: filterMappings (+ TTAA if configured) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  sampleDir="$(dirname "$pathPrefix")"
  echo "-- $name (in $sampleDir) --"
  if [ -n "$TTAA_BED" ]; then
    (cd "$sampleDir" && Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$sampleDir/" "$name" ttaaFilter "$TTAA_BED")
  else
    (cd "$sampleDir" && Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$sampleDir/" "$name")
  fi
done < "$PARSED/FASTQ_MAPS.tsv"
echo

# =============================================================================
# Stage 05: filterContaminants -- one run per contamination batch, from $WORK (reads each
# sample's .filtered_barcodes.bed via single_combined/'s symlinks back to run3//run5/).
# =============================================================================
cd "$WORK"
echo "=== Stage 05: filterContaminants ==="
BATCHES=$(awk -F'\t' '{print ($6==""?"default":$6)}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
while read -r batch; do
  [ -z "$batch" ] && continue
  FILELIST="$WORK/.pipeline_fileList.${batch}.txt"
  awk -F'\t' -v b="$batch" '($6==""?"default":$6)==b {print $1"\t"$5}' "$PARSED/FASTQ_MAPS.tsv" > "$FILELIST"
  echo "-- batch: $batch ($(wc -l < "$FILELIST") sample(s)) --"
  Rscript "$SCRIPTS/05_filterContaminants_uniform.R" "$WORK/" "$(basename "$FILELIST")" "$CONT_FACTOR"
done <<< "$BATCHES"
echo

# =============================================================================
# Stage 06: get_genuine_barcodes -- one run per distinct mapName
# =============================================================================
echo "=== Stage 06: get_genuine_barcodes ==="
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
while read -r mapName; do
  [ -z "$mapName" ] && continue
  manifest="$WORK/${mapName}.txt"
  awk -F'\t' -v m="$mapName" '$5==m{print $1".retained_barcodes.bed"}' "$PARSED/FASTQ_MAPS.tsv" > "$manifest"
  echo "-- $mapName ($(wc -l < "$manifest") sample(s)) --"
  Rscript "$SCRIPTS/06_get_genuine_barcodes.R" "$WORK/" "$(basename "$manifest")" "$BARCODE_LEN" "$DPNII_BED"
done <<< "$MAPNAMES"
echo

# =============================================================================
# Stage 07: normalizeTRIP -- one run per COMBOS row (combined_counts.txt)
# =============================================================================
echo "=== Stage 07: normalizeTRIP (combined_counts.txt) ==="
while IFS=$'\t' read -r mapName normSample expSample effLength revComp oneMismatch; do
  [ -z "$mapName" ] && continue
  echo "-- map=$mapName norm=$normSample exp=$expSample --"
  args=("$WORK/" "$WORK/" "${mapName}.txt" "$normSample" "$expSample" "$effLength" "$revComp")
  if [ "${oneMismatch:-}" = "yes" ]; then args+=("oneMismatch"); fi
  Rscript "$SCRIPTS/07_normalizeTRIP_MM.R" "${args[@]}"
done < "$PARSED/COMBOS.tsv"
echo

if [ -s "$PARSED/SPLICED_UNSPLICED.tsv" ]; then
  echo "=== Stage 14: combineSplicedUnsplicedPairs ==="
  cp "$PARSED/SPLICED_UNSPLICED.tsv" "$WORK/spliced_unspliced.txt"
  Rscript "$SCRIPTS/14_combineSplicedUnsplicedPairs_uniform.R" "$WORK/"
  echo
fi

# =============================================================================
# Stages 08-10: per-construct reference annotation
# =============================================================================
NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

echo "=== Stage 08: make_tracks (MANE/gencode/repeats intersects) ==="
bash "$SCRIPTS/08_make_tracks_genuine_barcodes.sh" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$REPEATS_BED"
echo

echo "=== Stage 09: gc_content (one run per construct) ==="
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "-- $mapName --"
  bash "$SCRIPTS/09_gc_content_single_genuine_barcodes.sh" "${mapName}.txt" "$CHROM_SIZES" "$GENOME_FASTA"
done <<< "$MAPNAMES"
echo

echo "=== Stage 10: table_annotation (annotated_tables.txt) ==="
Rscript "$SCRIPTS/10_table_annotation_genuine_barcodes.R" "$WORK/" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$TRANS2GENE" "$MANE_EXPRESSION" "$CHROMHMM_BED" "$SPIN_BED"
echo

echo "=== done ==="
