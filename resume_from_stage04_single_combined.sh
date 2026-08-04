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

# -----------------------------------------------------------------------------
# Job submission helpers -- same contract as run_pipeline.sh's submitOrRun/waitForJobs
# (one bsub per independent item, then block on the explicit job IDs). Stage 04's samples
# are fully independent -- each reads and writes only its own <name>.* files inside its own
# run directory -- so they are submitted concurrently instead of looped sequentially, which
# is what run_pipeline.sh already does for its own stage 04. Every later stage is submitted
# the same way -- one job per contamination batch (05), per construct (06, 09) and per COMBOS
# row (07) -- with a waitForJobs barrier between stages. The barriers are the real
# dependencies: 05 builds its cross-group contamination reference from EVERY sample's
# .filtered_barcodes.bed, 06 consumes every sample's .retained_barcodes.bed, and 09-11 need
# all of 07's combined_counts.txt. Within a stage the items are independent.
#
# waitForJobs polls each job by its own ID and requires DONE. Deliberately not
# `bsub -w "done(<jobname>)"`: with several jobs sharing a name LSF applies that to the most
# recently submitted one only, so a dependent stage can start while siblings are still
# writing (this silently produced a merge over stale inputs earlier in this project).
# -----------------------------------------------------------------------------
USE_LSF=${USE_LSF:-1}
LSF_QUEUE=${LSF_QUEUE:-gsla-cpu}
STAGE04_MEM=${STAGE04_MEM:-8000}
R_MOD=${R_MOD:-R/4.3.1-foss-2022b}
BEDTOOLS_MOD=${BEDTOOLS_MOD:-BEDTools/2.31.1-GCC-13.2.0}
JOBIDS=()

submitOrRun() {
  local logfile=$1 mem=$2 ncores=$3 modules=$4
  shift 4
  local preamble=""
  if [ "$modules" != "-" ] && [ -n "$modules" ]; then
    preamble="source /etc/profile.d/lmod.sh; module load $modules; "
  fi
  if [ "$USE_LSF" = "1" ]; then
    local out jobid
    out=$(bsub -q "$LSF_QUEUE" -R "rusage[mem=${mem}]" -n "$ncores" -o "$logfile" "${preamble}$*" 2>&1)
    jobid=$(echo "$out" | grep -oE '<[0-9]+>' | head -1 | tr -d '<>')
    if [ -z "$jobid" ]; then
      echo "ERROR: bsub submission failed:" >&2; echo "$out" >&2; exit 1
    fi
    JOBIDS+=("$jobid")
  else
    eval "${preamble}$*" > "$logfile" 2>&1
  fi
}

waitForJobs() {
  [ "$USE_LSF" = "1" ] || { JOBIDS=(); return 0; }
  [ ${#JOBIDS[@]} -eq 0 ] && return 0
  local failed=0 jid st
  for jid in "${JOBIDS[@]}"; do
    while true; do
      st=$(bjobs -noheader -o stat "$jid" 2>/dev/null || true)
      [ "$st" = "RUN" ] || [ "$st" = "PEND" ] || break
      sleep 15
    done
    st=$(bjobs -noheader -o stat "$jid" 2>/dev/null || true)
    if [ "$st" != "DONE" ]; then
      echo "ERROR: job $jid finished with status '${st:-UNKNOWN}' (expected DONE)" >&2
      failed=1
    fi
  done
  JOBIDS=()
  if [ "$failed" = "1" ]; then
    echo "One or more jobs in this stage failed -- aborting pipeline." >&2
    exit 1
  fi
}

# =============================================================================
# Stage 04: filterMappings (+ TTAA if configured) -- one job per sample, run in that
# sample's own run directory (dirname of its FASTQ_MAPS path_prefix), not in $WORK.
# =============================================================================
echo "=== Stage 04: filterMappings (+ TTAA if configured; one job per sample) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  sampleDir="$(dirname "$pathPrefix")"
  echo "-- $name (in $sampleDir) --"
  # Per-sample log: concurrent jobs would otherwise interleave into one stream.
  if [ -n "$TTAA_BED" ]; then
    submitOrRun "$sampleDir/${name}_04_filterMappings.log" "$STAGE04_MEM" 1 "$R_MOD" \
      "cd '$sampleDir' && Rscript '$SCRIPTS/04_filterMappings_TTAA.R' '$sampleDir/' '$name' ttaaFilter '$TTAA_BED'"
  else
    submitOrRun "$sampleDir/${name}_04_filterMappings.log" "$STAGE04_MEM" 1 "$R_MOD" \
      "cd '$sampleDir' && Rscript '$SCRIPTS/04_filterMappings_TTAA.R' '$sampleDir/' '$name'"
  fi
done < "$PARSED/FASTQ_MAPS.tsv"
waitForJobs
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
  submitOrRun "stage05_${batch}.log" 8000 1 "$R_MOD" \
    "cd '$WORK' && Rscript '$SCRIPTS/05_filterContaminants_uniform.R' '$WORK/' '$(basename "$FILELIST")' '$CONT_FACTOR'"
done <<< "$BATCHES"
waitForJobs
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
  submitOrRun "stage06_${mapName}.log" 8000 1 "$R_MOD" \
    "cd '$WORK' && Rscript '$SCRIPTS/06_get_genuine_barcodes.R' '$WORK/' '$(basename "$manifest")' '$BARCODE_LEN' '$DPNII_BED'"
done <<< "$MAPNAMES"
waitForJobs
echo

# =============================================================================
# Stage 07: normalizeTRIP -- one run per COMBOS row (combined_counts.txt)
# =============================================================================
echo "=== Stage 07: normalizeTRIP (combined_counts.txt) ==="
while IFS=$'\t' read -r mapName normSample expSample effLength revComp oneMismatch; do
  [ -z "$mapName" ] && continue
  echo "-- map=$mapName norm=$normSample exp=$expSample --"
  argstr="'$WORK/' '$WORK/' '${mapName}.txt' '$normSample' '$expSample' '$effLength' '$revComp'"
  if [ "${oneMismatch:-}" = "yes" ]; then argstr="$argstr 'oneMismatch'"; fi
  submitOrRun "stage07_${expSample}.log" 8000 1 "$R_MOD" \
    "cd '$WORK' && Rscript '$SCRIPTS/07_normalizeTRIP_MM.R' $argstr"
done < "$PARSED/COMBOS.tsv"
waitForJobs
echo

if [ -s "$PARSED/SPLICED_UNSPLICED.tsv" ]; then
  echo "=== Stage 08: combineSplicedUnsplicedPairs ==="
  cp "$PARSED/SPLICED_UNSPLICED.tsv" "$WORK/spliced_unspliced.txt"
  submitOrRun "stage08_combineSplicedUnspliced.log" 8000 1 "$R_MOD" \
    "cd '$WORK' && Rscript '$SCRIPTS/08_combineSplicedUnsplicedPairs_uniform.R' '$WORK/'"
  waitForJobs
  echo
fi

# =============================================================================
# Stages 09-11: per-construct reference annotation
# =============================================================================
NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

echo "=== Stage 09: make_tracks (MANE/gencode/repeats intersects) ==="
submitOrRun "stage09_make_tracks.log" 4000 1 "$BEDTOOLS_MOD" \
  "cd '$WORK' && bash '$SCRIPTS/09_make_tracks_genuine_barcodes.sh' '$(basename "$NAMELIST")' '$MANE_BED' '$MANE_INTRONS_BED' '$GENCODE_NONCODING_BED' '$REPEATS_BED'"
waitForJobs
echo

echo "=== Stage 10: gc_content (one run per construct) ==="
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "-- $mapName --"
  submitOrRun "stage10_${mapName}.log" 4000 1 "$BEDTOOLS_MOD" \
    "cd '$WORK' && bash '$SCRIPTS/10_gc_content_single_genuine_barcodes.sh' '${mapName}.txt' '$CHROM_SIZES' '$GENOME_FASTA'"
done <<< "$MAPNAMES"
waitForJobs
echo

echo "=== Stage 11: table_annotation (annotated_tables.txt + annotated_table_used.txt) ==="
# The 10th arg (norm manifest: mapName<TAB>normSample per COMBOS row) is what makes stage 11
# also emit <construct>.annotated_table_used.txt -- the same table restricted to barcodes with
# >=2 normalization reads. It was missing here, so this script regenerated .annotated_tables.txt
# but silently left .annotated_table_used.txt stale from whenever run_pipeline.sh last ran --
# a real hazard, since computeIntegrationStats.R and comparePairwiseSharedPositions_uniform.R
# read the _used files, which would then disagree with their own .annotated_tables.txt.
NORMMANIFEST="$WORK/.pipeline_norm_manifest.txt"
awk -F'\t' '{print $1"\t"$2}' "$PARSED/COMBOS.tsv" > "$NORMMANIFEST"
submitOrRun "stage11_table_annotation.log" 16000 1 "$R_MOD" \
  "cd '$WORK' && Rscript '$SCRIPTS/11_table_annotation_genuine_barcodes.R' '$WORK/' '$(basename "$NAMELIST")' '$MANE_BED' '$MANE_INTRONS_BED' '$GENCODE_NONCODING_BED' '$TRANS2GENE' '$MANE_EXPRESSION' '$CHROMHMM_BED' '$SPIN_BED' '$(basename "$NORMMANIFEST")'"
waitForJobs
echo

echo "=== done ==="
