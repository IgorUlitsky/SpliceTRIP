#!/bin/bash
# =============================================================================
# TRIP unified pipeline -- stages 12-13 only (metagene bigWig overlap + correlation report)
# =============================================================================
# Runs ONLY stage 12 (meta_bigwig_overlap) then stage 13 (plotMetageneCorrelations) against
# an already-completed dataset (stages 01-11 must have already produced that dataset's
# <mapName>.txt.retained_barcodes_genuine.bed maps and <expName>.combined_counts.txt files).
# run_pipeline.sh has no "start from stage N" mode, and re-running it with target=annotated
# on an already-processed dataset would needlessly redo cutadapt/bowtie2/etc for every
# sample. Sibling of run_pipeline.sh, which it leaves untouched -- shares its config
# format and its submitOrRun/waitForJobs helpers and stage 12/13 argument construction
# verbatim (copied rather than sourced, so this stays a single self-contained script like
# run_pipeline.sh itself).
#
# Usage
#   bash run_stages_12_13.sh <config.txt>
#
# Run from the dataset's working directory (same convention as run_pipeline.sh), e.g.:
#   cd /path/to/single_combined && bash <repo>/run_stages_12_13.sh <repo>/config/single_combined_config.txt
#
# Requires histone_marks.txt and speckle_names.txt (real, curated per-dataset files) in the
# working directory -- exits without running anything if either is missing, same gate
# run_pipeline.sh applies before its own stages 12-13.
# =============================================================================
set -euo pipefail

CONFIG=${1:?"usage: run_stages_12_13.sh <config.txt>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"
WORK=$(pwd)

echo "=== TRIP unified pipeline -- stages 12-13 only ==="
echo "config: $CONFIG"
echo "workdir: $WORK"
echo

if [ ! -f "$WORK/histone_marks.txt" ] || [ ! -f "$WORK/speckle_names.txt" ]; then
  echo "ERROR: histone_marks.txt and/or speckle_names.txt not found in $WORK -- stages 12-13 need both." >&2
  exit 1
fi

# -----------------------------------------------------------------------------
# Parse config into section files under a scratch dir, and source ### GLOBAL as env vars
# (same mechanism as run_pipeline.sh; only the sections stages 12-13 actually need).
# -----------------------------------------------------------------------------
PARSED=$(mktemp -d)
trap 'rm -rf "$PARSED"' EXIT

for section in GLOBAL FASTQ_MAPS COMBOS; do
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

: "${CHROM_SIZES:?CHROM_SIZES missing from ### GLOBAL}"
: "${SPECKLE_DIR:?SPECKLE_DIR missing from ### GLOBAL}"
: "${HISTONE_BW_DIR:?HISTONE_BW_DIR missing from ### GLOBAL}"
: "${HG38_TO_HG19_CHAIN:?HG38_TO_HG19_CHAIN missing from ### GLOBAL}"

echo "--- GLOBAL settings ---"
echo "SPECKLE_DIR=$SPECKLE_DIR  HISTONE_BW_DIR=$HISTONE_BW_DIR"
echo "CHROM_SIZES=$CHROM_SIZES  HG38_TO_HG19_CHAIN=$HG38_TO_HG19_CHAIN"
echo

# -----------------------------------------------------------------------------
# LSF parallel-submission helpers -- copied verbatim from run_pipeline.sh (see that
# script's own header comments for the reasoning behind waitForJobs' log-based, not
# bjobs-status-based, verdict).
# -----------------------------------------------------------------------------
USE_LSF=${USE_LSF:-1}
LSF_QUEUE=${LSF_QUEUE:-gsla-cpu}
declare -a JOBIDS=()
declare -A JOBLOG=()

R_MOD="R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1"
BEDTOOLS_MOD="BEDTools/2.31.1-GCC-13.3.0"
KENTUTILS_MOD="kentUtils/v377"

submitOrRun() {
  local logfile=$1 mem=$2 ncores=$3 modules=$4
  shift 4
  local preamble=""
  if [ "$modules" != "-" ] && [ -n "$modules" ]; then
    preamble="source /etc/profile.d/lmod.sh; module load $modules; "
  fi
  if [ "$USE_LSF" = "1" ]; then
    local out jobid
    rm -f "$logfile"
    out=$(bsub -q "$LSF_QUEUE" -R "select[type==X86_64 && model==AMD_EPYC] rusage[mem=${mem}]" -n "$ncores" -o "$logfile" "${preamble}$*" 2>&1)
    jobid=$(echo "$out" | grep -oE '<[0-9]+>' | head -1 | tr -d '<>')
    if [ -z "$jobid" ]; then
      echo "ERROR: bsub submission failed:" >&2
      echo "$out" >&2
      exit 1
    fi
    JOBIDS+=("$jobid")
    JOBLOG["$jobid"]="$logfile"
  else
    eval "${preamble}$*" > "$logfile" 2>&1
  fi
}

waitForJobs() {
  [ "$USE_LSF" = "1" ] || { JOBIDS=(); return 0; }
  [ ${#JOBIDS[@]} -eq 0 ] && return 0
  local failed=0 jid log active gone
  local GRACE_POLLS=8
  declare -A graceLeft=()
  for jid in "${JOBIDS[@]}"; do graceLeft["$jid"]=$GRACE_POLLS; done

  local -a pending=("${JOBIDS[@]}")
  while [ ${#pending[@]} -gt 0 ]; do
    active=$(bjobs -w 2>/dev/null | awk 'NR>1{print $1}' || true)
    local -a stillPending=()
    for jid in "${pending[@]}"; do
      log=${JOBLOG[$jid]:-}
      if grep -q "^Successfully completed\." "$log" 2>/dev/null; then
        continue
      fi
      if grep -qE "^(Exited with exit code|TERM_)" "$log" 2>/dev/null; then
        echo "ERROR: job $jid failed -- see $log" >&2
        grep -E "^(Exited with exit code|TERM_)" "$log" | head -2 >&2
        failed=1
        continue
      fi
      if grep -qx "$jid" <<< "$active"; then
        graceLeft["$jid"]=$GRACE_POLLS
        stillPending+=("$jid")
      else
        gone=$(( ${graceLeft[$jid]} - 1 ))
        graceLeft["$jid"]=$gone
        if [ "$gone" -le 0 ]; then
          echo "ERROR: job $jid is no longer active but wrote no completion marker to ${log:-<no log>}" >&2
          failed=1
        else
          stillPending+=("$jid")
        fi
      fi
    done
    pending=("${stillPending[@]}")
    [ ${#pending[@]} -gt 0 ] && sleep 15
  done

  JOBIDS=(); JOBLOG=()
  if [ "$failed" = "1" ]; then
    echo "One or more jobs in this stage failed -- aborting." >&2
    exit 1
  fi
}

# -----------------------------------------------------------------------------
# Stages 12-13 -- argument construction copied verbatim from run_pipeline.sh.
# -----------------------------------------------------------------------------
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
EXPLIST="$WORK/.pipeline_exp_list.txt"
awk -F'\t' '{print $1"\t"$3}' "$PARSED/COMBOS.tsv" > "$EXPLIST"
ALLMAPNAMES=$(echo "$MAPNAMES" | paste -sd,)
MAPTYPES="$WORK/.pipeline_map_types.txt"
printf '%s\n' "$MAPNAMES" > "$MAPTYPES"

echo "=== Stage 12: meta_bigwig_overlap (bedtools slop + bigWigAverageOverBed, one job per construct) ==="
# 12_meta_bigwig_overlap.sh's own header flags this: it's a straight sequential port of the
# original csh scripts, looping over every (construct x histone/speckle track) pair as one
# bigWigAverageOverBed call after another inside a single job -- 4 constructs x 59 tracks
# (single_combined) is 236 sequential external-tool calls in one job. The script itself
# already accepts a MAPTYPES_FILE with just one mapName, so no change to its own logic is
# needed -- fan out one job per construct here instead of a single job over the whole list,
# matching stage 10 (gc_content)'s existing one-job-per-construct convention elsewhere in
# this same file.
while read -r mapName; do
  [ -z "$mapName" ] && continue
  oneMapType="$WORK/.pipeline_map_types.${mapName}.txt"
  printf '%s\n' "$mapName" > "$oneMapType"
  submitOrRun "stage12_meta_bigwig_${mapName}.log" 8000 1 "$BEDTOOLS_MOD $KENTUTILS_MOD" bash "$SCRIPTS/12_meta_bigwig_overlap.sh" "$(basename "$oneMapType")" \
    "$HISTONE_BW_DIR" "$SPECKLE_DIR" "$CHROM_SIZES" "$HG38_TO_HG19_CHAIN"
done <<< "$MAPNAMES"
waitForJobs
echo

echo "=== Stage 13: plotMetageneCorrelations ==="
submitOrRun "stage13_metagene.log" 8000 1 "$R_MOD" Rscript "$SCRIPTS/13_plotMetageneCorrelations_uniform.R" "$WORK/" pearson "$ALLMAPNAMES" \
  "$(basename "$EXPLIST")" metagene_correlations
waitForJobs
echo

echo "=== done: metagene_correlations.{meta_corrs,summary_meta_corrs}.pdf written in $WORK/plots/ ==="
