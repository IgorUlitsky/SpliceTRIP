#!/bin/bash
# =============================================================================
# Backup every file that resume_from_stage04_introntrip_combined.sh will overwrite
# =============================================================================
# Usage
#   bash backup_before_resume_introntrip_combined.sh <config.txt> <work_dir> [backup_root]
#
# Scoped exactly to the outputs of stages 04-11 (nothing else is touched by the
# rerun, so nothing else is copied -- a full copy of IntronTRIP_combined/ would be 13 GB
# of mostly untouched cutadapt counts; this set is ~5 GB):
#
#   per ### FASTQ_MAPS sample, in that sample's OWN run directory (IntronTRIP_run2/3/4 --
#   this is where stages 04/05 actually write, IntronTRIP_combined/ only symlinks to them):
#     .filtered_barcodes.bed  .filtered_barcodes_poses.txt          (stage 04)
#     .retained_barcodes.bed  .retained_barcodes_counts.txt         (stage 05)
#
#   in <work_dir> (IntronTRIP_combined/), per construct <mapName>.txt:
#     .retained_barcodes_genuine.bed                                (stage 06)
#     .MANE_all/.MANE_exonic/.MANE_exonic_antisense/.MANE_all_antisense/
#     .MANE_introns_antisense/.MANE_introns_sense/.gencode_noncoding_all/.repeats .txt
#                                                                   (stage 10)
#     _tmp_extended.bed  _sequences.fasta  _gc_content.txt  _with_gc.bed
#                                                                   (stage 10)
#     .annotated_tables.txt  .annotated_table_used.txt              (stage 11)
#
#   in <work_dir>, per ### COMBOS expSample:  .bed  .combined_counts.txt   (stage 07)
#   in <work_dir>, per ### SPLICED_UNSPLICED outName:                      (stage 08)
#     .spliced_unspliced / .spliced_to_total / .spliced_and_unspliced_to_total
#       .combined_counts.txt
#
# Copies with `cp -L` (dereferencing) so a backed-up symlink becomes a real file -- the
# point of the backup is to survive the rerun overwriting the symlink's TARGET.
# Mirrors the layout of the existing backup_preContamFix_20260804_0728/ (one subdirectory
# per source directory).
# =============================================================================
set -euo pipefail

CONFIG=${1:?"usage: backup_before_resume_introntrip_combined.sh <config.txt> <work_dir> [backup_root]"}
WORK=${2:?"usage: backup_before_resume_introntrip_combined.sh <config.txt> <work_dir> [backup_root]"}
WORK=${WORK%/}
BACKUP_ROOT=${3:-"$(dirname "$WORK")/backup_IntronTRIP_preRerun_$(date +%Y%m%d_%H%M)"}

echo "=== Backup before IntronTRIP stage 04-11 rerun ==="
echo "config:      $CONFIG"
echo "work dir:    $WORK"
echo "backup root: $BACKUP_ROOT"
echo

PARSED=$(mktemp -d)
trap 'rm -rf "$PARSED"' EXIT
for section in FASTQ_MAPS COMBOS SPLICED_UNSPLICED; do
  awk -v s="### $section" '
    $0==s {grab=1; next}
    /^### /{grab=0}
    grab && NF && $0 !~ /^#/ {print}
  ' "$CONFIG" > "$PARSED/$section.tsv"
done

mkdir -p "$BACKUP_ROOT"

copied=0
missing=0
# copyOne <src> <dest_dir> -- skips (and reports) sources that don't exist yet, since a
# first-time stage output legitimately has nothing to back up.
copyOne() {
  local src=$1 dest=$2
  if [ -e "$src" ]; then
    cp -L -p "$src" "$dest/"
    copied=$((copied+1))
  else
    echo "   (absent, nothing to back up): $src"
    missing=$((missing+1))
  fi
}

# ---------------------------------------------------------------------------
# Stage 04/05 outputs, in each map sample's own run directory
# ---------------------------------------------------------------------------
echo "--- stage 04/05 per-sample outputs (run directories) ---"
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  sampleDir=$(dirname "$pathPrefix")
  destDir="$BACKUP_ROOT/$(basename "$sampleDir")"
  mkdir -p "$destDir"
  for suf in .filtered_barcodes.bed .filtered_barcodes_poses.txt \
             .retained_barcodes.bed .retained_barcodes_counts.txt; do
    copyOne "$sampleDir/$name$suf" "$destDir"
  done
done < "$PARSED/FASTQ_MAPS.tsv"

DEST_WORK="$BACKUP_ROOT/$(basename "$WORK")"
mkdir -p "$DEST_WORK"

# ---------------------------------------------------------------------------
# Stage 06/08/09/10 outputs, per construct, in the work dir
# ---------------------------------------------------------------------------
echo "--- stage 06/08/09/10 per-construct outputs ($(basename "$WORK")) ---"
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
while read -r mapName; do
  [ -z "$mapName" ] && continue
  for suf in .retained_barcodes_genuine.bed \
             .MANE_all.txt .MANE_exonic.txt .MANE_exonic_antisense.txt .MANE_all_antisense.txt \
             .MANE_introns_antisense.txt .MANE_introns_sense.txt .gencode_noncoding_all.txt .repeats.txt \
             _tmp_extended.bed _sequences.fasta _gc_content.txt _with_gc.bed \
             .annotated_tables.txt .annotated_table_used.txt; do
    copyOne "$WORK/${mapName}.txt$suf" "$DEST_WORK"
  done
done <<< "$MAPNAMES"

# ---------------------------------------------------------------------------
# Stage 07 outputs, per COMBOS expSample (deduplicated -- an expSample can appear in
# several COMBOS rows)
# ---------------------------------------------------------------------------
echo "--- stage 07 per-expSample outputs ($(basename "$WORK")) ---"
while read -r expSample; do
  [ -z "$expSample" ] && continue
  copyOne "$WORK/$expSample.bed" "$DEST_WORK"
  copyOne "$WORK/$expSample.combined_counts.txt" "$DEST_WORK"
done < <(awk -F'\t' '{print $3}' "$PARSED/COMBOS.tsv" | sort -u)

# ---------------------------------------------------------------------------
# Stage 08 outputs, per SPLICED_UNSPLICED row
# ---------------------------------------------------------------------------
echo "--- stage 08 splicing-ratio outputs ($(basename "$WORK")) ---"
while IFS=$'\t' read -r outName spliced unspliced total; do
  [ -z "$outName" ] && continue
  for suf in .spliced_unspliced.combined_counts.txt \
             .spliced_to_total.combined_counts.txt \
             .spliced_and_unspliced_to_total.combined_counts.txt; do
    copyOne "$WORK/$outName$suf" "$DEST_WORK"
  done
done < "$PARSED/SPLICED_UNSPLICED.tsv"

echo
echo "=== backup done: $copied file(s) copied, $missing absent ==="
du -sh "$BACKUP_ROOT"
echo "$BACKUP_ROOT"
