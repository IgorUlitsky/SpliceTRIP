#!/bin/bash
# =============================================================================
# One-off resume script: re-run IntronTRIP_combined from stage 04 onward
# (stages 04 -> 05 -> 06 -> 07 -> 14 -> 08 -> 09 -> 10), skipping stages 01-03 since their
# outputs (.barcodes.bed / .barcodes_poses.txt) already exist for all 32 ### FASTQ_MAPS
# samples. Same motivation as resume_from_stage04_single_combined.sh: the stage 04/05
# outputs currently on disk (IntronTRIP_run2/3/4/*.filtered_barcodes.bed dated 2026-07-24/25,
# *.retained_barcodes.bed dated 2026-07-30) predate both this repo's stage 04 barcode-dedup
# logic and the 2026-08-03 stage 05 fix ("key the contamination reference per barcode"), so
# everything from stage 04 down needs regenerating.
#
# Two structural differences from run_pipeline.sh, both matching how the REAL IntronTRIP
# data is laid out on disk:
#
#  1. Stage 03's outputs live in each sample's OWN run directory (IntronTRIP_run2/, _run3/,
#     _run4/), not in IntronTRIP_combined/ -- so stage 04 runs once per sample, cd'd into
#     that sample's run directory. Stages 05+ run from IntronTRIP_combined/, where the
#     cross-sample manifests and symlinks live. Same split as
#     resume_from_stage04_single_combined.sh.
#
#  2. IntronTRIP_combined/ is MISSING some of the per-sample symlinks that stage 05 reads
#     and writes through (verified 2026-08-04: no *.filtered_barcodes_poses.txt symlinks at
#     all, and the 12 IntronTRIP_run4 B1/C1/D1 samples additionally have no
#     *.filtered_barcodes.bed symlink) -- unlike single_combined/, which has the full set.
#     Without them stage 05 either errors on a missing input or writes a real file into
#     IntronTRIP_combined/ that diverges from the run directory's copy. linkStageFiles()
#     below creates any missing link before stage 05 runs, pointing at the run directory,
#     matching single_combined/'s layout exactly.
#
# Everything else reuses run_pipeline.sh's own config parsing and per-stage command blocks
# verbatim (including its LSF submitOrRun/waitForJobs parallelization, one job per
# sample/batch/construct/COMBOS row), so behavior stays identical to the real pipeline.
# USE_LSF=0 in the environment falls back to sequential local execution.
#
# Usage
#   bash resume_from_stage04_introntrip_combined.sh <config.txt> <work_dir> [threads]
# e.g.
#   bash resume_from_stage04_introntrip_combined.sh \
#     config/introntrip_combined_config.txt \
#     /home/projects/ulitsky/shared/data/human/TRIP/IntronTRIP_combined/
#
# Run backup_before_resume_introntrip_combined.sh FIRST -- every stage below overwrites its
# previous output in place.
# =============================================================================
set -euo pipefail

CONFIG=${1:?"usage: resume_from_stage04_introntrip_combined.sh <config.txt> <work_dir> [start_stage]"}
WORK=${2:?"usage: resume_from_stage04_introntrip_combined.sh <config.txt> <work_dir> [start_stage]"}
# start_stage (default 04) -- skip stages before it, for picking a run back up after an
# interruption. Only 04 and 05 are offered: everything from 05 on runs as one dependent
# chain from $WORK, whereas 04 is the one stage with a different working directory per
# sample, so "did stage 04 finish?" is the only checkpoint worth resuming from.
START_STAGE=${3:-04}
WORK=${WORK%/}
case "$START_STAGE" in
  04|4|05|5) ;;
  *) echo "unknown start_stage '$START_STAGE' (expected 04 or 05)" >&2; exit 1;;
esac
[ "$START_STAGE" = "4" ] && START_STAGE=04
[ "$START_STAGE" = "5" ] && START_STAGE=05
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"

echo "=== IntronTRIP_combined: resume from stage 04 (04-11) ==="
echo "config:   $CONFIG"
echo "work dir: $WORK"
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

echo "BARCODE_LEN=$BARCODE_LEN  CONT_FACTOR=$CONT_FACTOR  TTAA_BED=${TTAA_BED:-<off>}  DPNII_BED=${DPNII_BED:-<script default>}"
echo

# -----------------------------------------------------------------------------
# LSF helpers -- copied from run_pipeline.sh, plus the "-q gsla-cpu" this cluster's jobs
# must use (LSF's own configured default queue here is `short`, which is the wrong one).
# -----------------------------------------------------------------------------
USE_LSF=${USE_LSF:-1}
QUEUE=${QUEUE:-gsla-cpu}
declare -a JOBIDS=()
declare -A JOBLOG=()

R_MOD="R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1"
BEDTOOLS_MOD="BEDTools/2.31.1-GCC-13.3.0"

# submitOrRun <logfile> <mem_mb> <n_cores> <modules|-> <workdir> <command...>
# Unlike run_pipeline.sh's version this takes an explicit working directory, because stage
# 04 runs per sample inside that sample's own run directory while every later stage runs
# from $WORK.
submitOrRun() {
  local logfile=$1 mem=$2 ncores=$3 modules=$4 rundir=$5
  shift 5
  local preamble="cd $rundir; "
  if [ "$modules" != "-" ] && [ -n "$modules" ]; then
    preamble="source /etc/profile.d/lmod.sh; module load $modules; cd $rundir; "
  fi
  if [ "$USE_LSF" = "1" ]; then
    local out jobid
    # bsub -o APPENDS to an existing file (only -oo overwrites), so a stale log from an
    # earlier run would still carry its old "Successfully completed."/"Exited with exit
    # code" epilogue -- which is exactly what waitForJobs reads to decide this job's fate.
    # Remove it first so every marker found belongs to THIS submission.
    rm -f "$logfile"
    # Compiled-binary modules (bedtools) live under an AMD-specific easybuild tree not
    # visible on non-AMD_EPYC hosts, so pin every job to that host class.
    out=$(bsub -q "$QUEUE" -R "select[type==X86_64 && model==AMD_EPYC] rusage[mem=${mem}]" -n "$ncores" -o "$logfile" "${preamble}$*" 2>&1)
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

# waitForJobs -- blocks until every job in JOBIDS has ended, then verifies each one from the
# epilogue LSF writes into that job's own -o log ("Successfully completed." vs "Exited with
# exit code N"). Aborts the whole run if any job failed.
#
# Deliberately does NOT decide anything from `bjobs -o stat <id>`, which is what
# run_pipeline.sh's version does and what aborted the 2026-08-04 08:39 run after all 32
# stage-04 jobs had in fact succeeded: `bjobs <id>` yields an empty status both once LSF has
# purged the finished job's record and whenever mbatchd is transiently unreachable ("mbatchd
# is down or the ports are inaccessable", seen on this cluster at submission time the same
# morning) -- and empty was being mapped to UNKNOWN, i.e. failure. The job's own log is
# authoritative and never expires; `bjobs -w` (active listing, RUN/PEND only) is used just as
# the cheap "is it still going?" signal, per the same reasoning.
waitForJobs() {
  [ "$USE_LSF" = "1" ] || { JOBIDS=(); return 0; }
  [ ${#JOBIDS[@]} -eq 0 ] && return 0
  local failed=0 jid log active still gone

  # An ended job normally has its epilogue within seconds; allow GRACE_POLLS x 15 s of
  # "absent from bjobs but no marker yet" (slow shared-filesystem flush) before giving up on
  # it. A job that vanished without ever writing an epilogue is a genuine failure.
  local GRACE_POLLS=8
  declare -A graceLeft=()
  for jid in "${JOBIDS[@]}"; do graceLeft["$jid"]=$GRACE_POLLS; done

  local -a pending=("${JOBIDS[@]}")
  while [ ${#pending[@]} -gt 0 ]; do
    # `bjobs -w` lists only RUN/PEND jobs; it exits non-zero with an empty list when nothing
    # is active, which is indistinguishable from a transient failure -- hence the log check
    # below is what actually decides, and absence here only starts the grace countdown.
    active=$(bjobs -w 2>/dev/null | awk 'NR>1{print $1}' || true)
    local -a stillPending=()
    for jid in "${pending[@]}"; do
      log=${JOBLOG[$jid]:-}
      if grep -q "^Successfully completed\." "$log" 2>/dev/null; then
        continue                                   # done, verified
      fi
      if grep -qE "^(Exited with exit code|TERM_)" "$log" 2>/dev/null; then
        echo "ERROR: job $jid failed -- see $log" >&2
        grep -E "^(Exited with exit code|TERM_)" "$log" | head -2 >&2
        failed=1
        continue
      fi
      if grep -qx "$jid" <<< "$active"; then
        graceLeft["$jid"]=$GRACE_POLLS             # still active, reset the countdown
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

# =============================================================================
# Stage 04: filterMappings (+ TTAA if configured) -- one job per sample, each cd'd into
# that sample's own run directory (where stage 03's .barcodes.bed/.barcodes_poses.txt are).
# =============================================================================
if [ "$START_STAGE" = "05" ]; then
  echo "=== Stage 04 skipped (start_stage=05) ==="
  # Guard the skip: stage 05 silently reads whatever .filtered_barcodes* happens to be on
  # disk, so confirm every sample really does have a verified-complete stage 04 log rather
  # than trusting the caller's "it finished".
  incomplete=0
  while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
    [ -z "$name" ] && continue
    l="$(dirname "$pathPrefix")/${name}_04_filterMappings.log"
    if ! grep -q "^Successfully completed\." "$l" 2>/dev/null; then
      echo "   MISSING/incomplete stage 04 for $name ($l)" >&2
      incomplete=$((incomplete+1))
    fi
  done < "$PARSED/FASTQ_MAPS.tsv"
  if [ "$incomplete" -gt 0 ]; then
    echo "$incomplete sample(s) have no verified stage 04 output -- rerun with start_stage=04." >&2
    exit 1
  fi
  echo "   verified: all $(grep -c . "$PARSED/FASTQ_MAPS.tsv") sample(s) have a completed stage 04"
  echo
else
echo "=== Stage 04: filterMappings (+ TTAA if configured; one job per sample) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  sampleDir="$(dirname "$pathPrefix")"
  echo "-- $name (in $sampleDir) --"
  # 16 GB: the IntronTRIP_run4 B1/C1/D1 map libraries are ~25 M-row .barcodes.bed / ~39 M-row
  # .barcodes_poses.txt, an order of magnitude bigger than the run2/run3 ones.
  if [ -n "$TTAA_BED" ]; then
    submitOrRun "$sampleDir/${name}_04_filterMappings.log" 16000 1 "$R_MOD" "$sampleDir" \
      Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$sampleDir/" "$name" ttaaFilter "$TTAA_BED"
  else
    submitOrRun "$sampleDir/${name}_04_filterMappings.log" 16000 1 "$R_MOD" "$sampleDir" \
      Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$sampleDir/" "$name"
  fi
done < "$PARSED/FASTQ_MAPS.tsv"
waitForJobs
echo
fi

# =============================================================================
# Link each sample's stage 04/05 files into $WORK if they aren't already linked (see header
# note 2). Stage 05 reads .filtered_barcodes{.bed,_poses.txt} from $WORK and writes
# .retained_barcodes{.bed,_counts.txt} there; with these links in place both resolve to the
# sample's own run directory, exactly as single_combined/ is already set up.
# =============================================================================
echo "=== Linking per-sample stage 04/05 files into $WORK ==="
linked=0
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  sampleDir="$(cd "$(dirname "$pathPrefix")" && pwd)"
  for suf in .filtered_barcodes.bed .filtered_barcodes_poses.txt \
             .retained_barcodes.bed .retained_barcodes_counts.txt; do
    # -e follows symlinks, so an existing-but-dangling link is repointed rather than kept.
    if [ ! -e "$WORK/$name$suf" ]; then
      ln -sfn "$sampleDir/$name$suf" "$WORK/$name$suf"
      echo "   linked $name$suf -> $sampleDir/"
      linked=$((linked+1))
    fi
  done
done < "$PARSED/FASTQ_MAPS.tsv"
echo "   $linked link(s) created"
echo

cd "$WORK"

# =============================================================================
# Stage 05: filterContaminants -- one job per contamination batch (run2/run3/run4)
# =============================================================================
echo "=== Stage 05: filterContaminants (one job per batch) ==="
BATCHES=$(awk -F'\t' '{print ($6==""?"default":$6)}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
while read -r batch; do
  [ -z "$batch" ] && continue
  FILELIST="$WORK/.pipeline_fileList.${batch}.txt"
  awk -F'\t' -v b="$batch" '($6==""?"default":$6)==b {print $1"\t"$5}' "$PARSED/FASTQ_MAPS.tsv" > "$FILELIST"
  echo "-- batch: $batch ($(wc -l < "$FILELIST") sample(s)) --"
  # 32 GB: this stage holds every sample in the batch in memory at once, and the run4 batch
  # is 12 samples x ~48 M-row .filtered_barcodes_poses.txt.
  submitOrRun "$WORK/stage05_${batch}.log" 32000 1 "$R_MOD" "$WORK" \
    Rscript "$SCRIPTS/05_filterContaminants_uniform.R" "$WORK/" "$(basename "$FILELIST")" "$CONT_FACTOR"
done <<< "$BATCHES"
waitForJobs
echo

# =============================================================================
# Stage 06: get_genuine_barcodes -- one job per distinct mapName (construct)
# =============================================================================
echo "=== Stage 06: get_genuine_barcodes (one job per construct) ==="
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
while read -r mapName; do
  [ -z "$mapName" ] && continue
  manifest="$WORK/${mapName}.txt"
  awk -F'\t' -v m="$mapName" '$5==m{print $1".retained_barcodes.bed"}' "$PARSED/FASTQ_MAPS.tsv" > "$manifest"
  echo "-- $mapName ($(wc -l < "$manifest") sample(s)) --"
  submitOrRun "$WORK/stage06_${mapName}.log" 16000 1 "$R_MOD" "$WORK" \
    Rscript "$SCRIPTS/06_get_genuine_barcodes.R" "$WORK/" "$(basename "$manifest")" "$BARCODE_LEN" "$DPNII_BED"
done <<< "$MAPNAMES"
waitForJobs
echo

# =============================================================================
# Stage 07: normalizeTRIP -- one job per ### COMBOS row (48 rows)
# =============================================================================
echo "=== Stage 07: normalizeTRIP (combined_counts.txt, one job per COMBOS row) ==="
while IFS=$'\t' read -r mapName normSample expSample effLength revComp oneMismatch; do
  [ -z "$mapName" ] && continue
  echo "-- map=$mapName norm=$normSample exp=$expSample --"
  args=("$WORK/" "$WORK/" "${mapName}.txt" "$normSample" "$expSample" "$effLength" "$revComp")
  if [ "${oneMismatch:-}" = "yes" ]; then args+=("oneMismatch"); fi
  submitOrRun "$WORK/stage07_${expSample}.log" 8000 1 "$R_MOD" "$WORK" \
    Rscript "$SCRIPTS/07_normalizeTRIP_MM.R" "${args[@]}"
done < "$PARSED/COMBOS.tsv"
waitForJobs
echo

# =============================================================================
# Stage 08: combineSplicedUnsplicedPairs -- single job over all ### SPLICED_UNSPLICED rows,
# depends on stage 07's combined_counts.txt. Same position in the order as run_pipeline.sh
# (right after 07, before 09-11).
# =============================================================================
if [ -s "$PARSED/SPLICED_UNSPLICED.tsv" ]; then
  echo "=== Stage 08: combineSplicedUnsplicedPairs ($(wc -l < "$PARSED/SPLICED_UNSPLICED.tsv") triplet(s)) ==="
  cp "$PARSED/SPLICED_UNSPLICED.tsv" "$WORK/spliced_unspliced.txt"
  submitOrRun "$WORK/stage08_combineSplicedUnspliced.log" 8000 1 "$R_MOD" "$WORK" \
    Rscript "$SCRIPTS/08_combineSplicedUnsplicedPairs_uniform.R" "$WORK/"
  waitForJobs
  echo
fi

# =============================================================================
# Stages 09-11: per-construct reference annotation (only needs stage 06's output)
# =============================================================================
NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

echo "=== Stage 09: make_tracks (MANE/gencode/repeats intersects) ==="
submitOrRun "$WORK/stage09_make_tracks.log" 8000 1 "$BEDTOOLS_MOD" "$WORK" \
  bash "$SCRIPTS/09_make_tracks_genuine_barcodes.sh" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$REPEATS_BED"
waitForJobs
echo

echo "=== Stage 10: gc_content (one job per construct) ==="
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "-- $mapName --"
  submitOrRun "$WORK/stage10_${mapName}.log" 8000 1 "$BEDTOOLS_MOD" "$WORK" \
    bash "$SCRIPTS/10_gc_content_single_genuine_barcodes.sh" "${mapName}.txt" "$CHROM_SIZES" "$GENOME_FASTA"
done <<< "$MAPNAMES"
waitForJobs
echo

echo "=== Stage 11: table_annotation (annotated_tables.txt + annotated_table_used.txt) ==="
NORMMANIFEST="$WORK/.pipeline_norm_manifest.txt"
awk -F'\t' '{print $1"\t"$2}' "$PARSED/COMBOS.tsv" > "$NORMMANIFEST"
submitOrRun "$WORK/stage11_table_annotation.log" 32000 1 "$R_MOD" "$WORK" \
  Rscript "$SCRIPTS/11_table_annotation_genuine_barcodes.R" "$WORK/" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$TRANS2GENE" "$MANE_EXPRESSION" \
  "$CHROMHMM_BED" "$SPIN_BED" "$(basename "$NORMMANIFEST")"
waitForJobs
echo

echo "=== done: stages 04-11 rerun for all $(echo "$MAPNAMES" | wc -l) constructs ==="
