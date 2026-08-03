#!/bin/bash
# One-off resume script: re-run single_combined from stage 07 onward, skipping 01-06.
# Needed after fixing 4 mislabeled rows in combos.txt (mapName=5ss_noP instead of CPA_noP
# for the 3intron_CPA_noP_* rows) -- stages 01-06 only depend on FASTQ_MAPS/barcode
# mapping, not COMBOS, so their (already-correct, freshly regenerated) outputs don't need
# redoing. Stages 08-10 also don't depend on COMBOS/combined_counts.txt (see
# run_pipeline.sh's own comment: "independent of stage 07, runs alongside it"), so this
# only needs to touch stage 07 (+ stage 14 if configured) -- included anyway for
# completeness/parity with a real target=annotated run, at negligible extra cost since 08-10
# don't depend on anything this changed.
set -euo pipefail

CONFIG=${1:?"usage: resume_from_stage07_single_combined.sh <config.txt> <single_combined_dir>"}
WORK=${2:?"usage: resume_from_stage07_single_combined.sh <config.txt> <single_combined_dir>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"

echo "=== Resume from stage 07 ==="
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

cd "$WORK"

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

echo "=== done ==="
