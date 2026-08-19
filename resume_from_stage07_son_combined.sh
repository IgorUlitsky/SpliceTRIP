#!/bin/bash
# One-off resume script: re-run SON_combined from stage 07 onward, skipping 01-06 and
# skipping the optional metagene stages 12-13.
#
# WHY (2026-08-19): SON_combined's SON_A1/SON_E1/SON_B3_9 genuine-barcode maps were rebuilt
# with the current stage 05, whose 2026-08-03 fix (commit 4ae99cd, per-barcode contamination
# reference) drops 3.0-7.5% of previously-retained barcodes per construct -- the old maps
# came from the pre-fix standalone filterContaminants_uniform.R. Everything downstream of
# stage 06 was therefore computed against stale maps and needs regenerating. SON_B1 was
# already produced post-fix and its stage 06 output is bit-identical, but its COMBOS rows are
# re-run here too: stage 08 loops over every ### SPLICED_UNSPLICED triplet internally and
# stage 11 processes every construct in one call, so a uniform rerun is both simpler and less
# error-prone than a partial one, and B1's rows recompute to identical content from unchanged
# inputs.
#
# Stages 09-11 do NOT depend on stage 07 (see run_pipeline.sh: "independent of stage 07"; in
# particular stage 11 reads <normSample>.counts.txt.gz directly, not combined_counts.txt), so
# the 09->10->11 chain is submitted up front with LSF -w dependencies and runs concurrently
# with the stage 07 fan-out instead of queueing behind it.
#
# Stages 12-13 (metagene bigWig overlap + its correlation report) are deliberately NOT run.
set -euo pipefail

CONFIG=${1:?"usage: resume_from_stage07_son_combined.sh <config.txt> <son_combined_dir>"}
WORK=${2:?"usage: resume_from_stage07_son_combined.sh <config.txt> <son_combined_dir>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"

USE_LSF=${USE_LSF:-1}
LSF_QUEUE=${LSF_QUEUE:-gsla-cpu}
R_MOD="R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1"
BEDTOOLS_MOD="BEDTools/2.31.1-GCC-13.3.0"

declare -a JOBIDS=()
declare -A JOBLOG=()

# submitOrRun <logfile> <mem_mb> <n_cores> <modules|-> [-w <lsf_dep>] <command...>
# Same contract as run_pipeline.sh's submitOrRun (AMD_EPYC pinning, module preamble, log-
# verified completion), plus an optional -w LSF dependency so independent chains can be
# submitted up front rather than serialized behind a bash-level barrier.
submitOrRun() {
  local logfile=$1 mem=$2 ncores=$3 modules=$4
  shift 4
  local dep=""
  if [ "${1:-}" = "-w" ]; then dep=$2; shift 2; fi
  local preamble=""
  if [ "$modules" != "-" ] && [ -n "$modules" ]; then
    preamble="source /etc/profile.d/lmod.sh; module load $modules; "
  fi
  if [ "$USE_LSF" = "1" ]; then
    local out jobid
    rm -f "$logfile"
    local -a depArgs=()
    [ -n "$dep" ] && depArgs=(-w "$dep")
    out=$(bsub -q "$LSF_QUEUE" -R "select[type==X86_64 && model==AMD_EPYC] rusage[mem=${mem}]" \
          -n "$ncores" -o "$logfile" "${depArgs[@]}" "${preamble}$*" 2>&1)
    jobid=$(echo "$out" | grep -oE '<[0-9]+>' | head -1 | tr -d '<>')
    if [ -z "$jobid" ]; then
      echo "ERROR: bsub submission failed:" >&2
      echo "$out" >&2
      exit 1
    fi
    JOBIDS+=("$jobid")
    JOBLOG[$jobid]="$(cd "$(dirname "$logfile")" && pwd)/$(basename "$logfile")"
    LAST_JOBID=$jobid
  else
    eval "${preamble}$*"
    LAST_JOBID=""
  fi
}

# waitForJobs [jobid...] -- verbatim behavior of run_pipeline.sh's waitForJobs (each job's
# fate read from its own LSF log epilogue, not `bjobs -o stat`), except that an explicit
# list of job ids can be passed to wait on a subset instead of everything submitted so far.
waitForJobs() {
  [ "$USE_LSF" = "1" ] || { JOBIDS=(); return 0; }
  local -a target=()
  if [ $# -gt 0 ]; then target=("$@"); else target=("${JOBIDS[@]:-}"); fi
  [ ${#target[@]} -eq 0 ] && return 0
  [ -z "${target[0]}" ] && return 0
  local failed=0 jid log active gone
  local GRACE_POLLS=8
  declare -A graceLeft=()
  for jid in "${target[@]}"; do graceLeft["$jid"]=$GRACE_POLLS; done
  local -a pending=("${target[@]}")
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
  if [ $# -eq 0 ]; then JOBIDS=(); fi
  if [ "$failed" = "1" ]; then
    echo "One or more jobs in this stage failed -- aborting." >&2
    exit 1
  fi
}

echo "=== Resume SON_combined from stage 07 (stages 12-13 skipped) ==="
echo "config: $CONFIG"
echo "workdir: $WORK"
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

cd "$WORK"

# -----------------------------------------------------------------------------
# Stages 09-11 first (submitted, not waited on): independent of stage 07, chained to each
# other by explicit LSF job id so they progress while stage 07's fan-out runs.
# -----------------------------------------------------------------------------
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

echo "=== Stage 09: make_tracks (MANE/gencode/repeats intersects) ==="
submitOrRun "stage09_make_tracks.log" 4000 1 "$BEDTOOLS_MOD" \
  bash "$SCRIPTS/09_make_tracks_genuine_barcodes.sh" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$REPEATS_BED"
JID09=$LAST_JOBID
echo "   job $JID09"

echo "=== Stage 10: gc_content (one job per construct, after stage 09) ==="
JIDS10=()
while read -r mapName; do
  [ -z "$mapName" ] && continue
  submitOrRun "stage10_${mapName}.log" 4000 1 "$BEDTOOLS_MOD" -w "ended($JID09)" \
    bash "$SCRIPTS/10_gc_content_single_genuine_barcodes.sh" "${mapName}.txt" "$CHROM_SIZES" "$GENOME_FASTA"
  JIDS10+=("$LAST_JOBID")
  echo "   $mapName -> job $LAST_JOBID"
done <<< "$MAPNAMES"

echo "=== Stage 11: table_annotation (after all of stage 10) ==="
NORMMANIFEST="$WORK/.pipeline_norm_manifest.txt"
awk -F'\t' '{print $1"\t"$2}' "$PARSED/COMBOS.tsv" > "$NORMMANIFEST"
DEP11=$(printf "ended(%s) && " "${JIDS10[@]}"); DEP11=${DEP11% && }
submitOrRun "stage11_table_annotation.log" 16000 1 "$R_MOD" -w "$DEP11" \
  Rscript "$SCRIPTS/11_table_annotation_genuine_barcodes.R" "$WORK/" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$TRANS2GENE" "$MANE_EXPRESSION" \
  "$CHROMHMM_BED" "$SPIN_BED" "$(basename "$NORMMANIFEST")"
JID11=$LAST_JOBID
echo "   job $JID11"
echo

# -----------------------------------------------------------------------------
# Stage 07: one job per ### COMBOS row, then stage 08 once they have all landed.
# -----------------------------------------------------------------------------
echo "=== Stage 07: normalizeTRIP (combined_counts.txt, one job per COMBOS row) ==="
JIDS07=()
while IFS=$'\t' read -r mapName normSample expSample effLength revComp oneMismatch; do
  [ -z "$mapName" ] && continue
  args=("$WORK/" "$WORK/" "${mapName}.txt" "$normSample" "$expSample" "$effLength" "$revComp")
  if [ "${oneMismatch:-}" = "yes" ]; then args+=("oneMismatch"); fi
  submitOrRun "stage07_${expSample}.log" 8000 1 "$R_MOD" Rscript "$SCRIPTS/07_normalizeTRIP_MM.R" "${args[@]}"
  JIDS07+=("$LAST_JOBID")
done < "$PARSED/COMBOS.tsv"
echo "   ${#JIDS07[@]} jobs submitted"
waitForJobs "${JIDS07[@]}"
echo "   stage 07 complete"
echo

if [ -s "$PARSED/SPLICED_UNSPLICED.tsv" ]; then
  echo "=== Stage 08: combineSplicedUnsplicedPairs ($(wc -l < "$PARSED/SPLICED_UNSPLICED.tsv") triplets) ==="
  cp "$PARSED/SPLICED_UNSPLICED.tsv" "$WORK/spliced_unspliced.txt"
  submitOrRun "stage08_combineSplicedUnspliced.log" 8000 1 "$R_MOD" \
    Rscript "$SCRIPTS/08_combineSplicedUnsplicedPairs_uniform.R" "$WORK/"
  JID08=$LAST_JOBID
  echo "   job $JID08"
  waitForJobs "$JID08"
  echo "   stage 08 complete"
  echo
fi

echo "=== waiting on the stage 09-11 chain ==="
waitForJobs "$JID09" "${JIDS10[@]}" "$JID11"
echo "   stages 09-11 complete"
echo
echo "=== done (stages 07-11; 12-13 intentionally skipped) ==="
