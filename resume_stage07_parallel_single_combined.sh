#!/bin/bash
# Stage 07 (normalizeTRIP -> combined_counts.txt) is embarrassingly parallel: each COMBOS
# row reads its own map/norm/exp inputs and writes its own <expSample>.combined_counts.txt,
# no shared state between rows. Submits one bsub job per row instead of running them in a
# sequential loop.
set -euo pipefail

CONFIG=${1:?"usage: resume_stage07_parallel_single_combined.sh <config.txt> <single_combined_dir>"}
WORK=${2:?"usage: resume_stage07_parallel_single_combined.sh <config.txt> <single_combined_dir>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"

PARSED=$(mktemp -d)
trap 'rm -rf "$PARSED"' EXIT
awk '
  $0=="### COMBOS" {grab=1; next}
  /^### /{grab=0}
  grab && NF && $0 !~ /^#/ {print}
' "$CONFIG" > "$PARSED/COMBOS.tsv"

cd "$WORK"
while IFS=$'\t' read -r mapName normSample expSample effLength revComp oneMismatch; do
  [ -z "$mapName" ] && continue
  args=("$WORK/" "$WORK/" "${mapName}.txt" "$normSample" "$expSample" "$effLength" "$revComp")
  if [ "${oneMismatch:-}" = "yes" ]; then args+=("oneMismatch"); fi
  echo "submitting: map=$mapName norm=$normSample exp=$expSample"
  bsub -q "${LSF_QUEUE:-gsla-cpu}" -R "rusage[mem=8000]" -o "${WORK}/stage07_${expSample}.log" \
    "source /etc/profile.d/lmod.sh; module load R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1; Rscript $SCRIPTS/07_normalizeTRIP_MM.R ${args[*]}"
done < "$PARSED/COMBOS.tsv"
