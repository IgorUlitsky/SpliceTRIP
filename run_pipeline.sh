#!/bin/bash
# =============================================================================
# TRIP unified pipeline — orchestrator
# =============================================================================
# Reads a config text file (see config/example_config.txt) and runs stages 01-13
# (cutadapt -> bowtie2 -> analyzeTRIPMapping -> filterMappings[+TTAA] ->
# filterContaminants -> get_genuine_barcodes -> normalizeTRIP[+oneMismatch] ->
# combineSplicedUnsplicedPairs[if configured] -> make_tracks -> gc_content ->
# table_annotation -> meta_bigwig_overlap[if configured] ->
# plotMetageneCorrelations[if configured]) in order. comparePairwiseSharedPositions /
# annotateTRIP_positions / annotateTRIPGeneBodies / gc_content_report /
# compareReplicatesAndMultiIntegrations are visualization/QC scripts, not pipeline stages --
# run them by hand against stage 07/11's output (see README's "Visualization scripts"
# section), not orchestrated here.
#
# Usage
#   bash run_pipeline.sh <config.txt> [target] [threads]
#
#   target (default: combined_counts) — how far to run:
#     cutadapt          stage 01 only
#     mapping           01-02 (Map-mode samples only)
#     genuine_barcodes  01-06
#     combined_counts   01-07 (everything needed for combined_counts.txt) ★ pipeline's main
#                        deliverable, plus stage 08 if ### SPLICED_UNSPLICED has rows
#     annotated         01-13 (adds per-construct annotated_tables.txt (+
#                        annotated_table_used.txt, the same table restricted to barcodes
#                        whose normalization reads SUMMED across that construct's ### COMBOS
#                        normSamples reach >=2 -- a pooled total, NOT ">=2 in at least one
#                        library"; see scripts/11_table_annotation_genuine_barcodes.R's
#                        header note) -- MANE/gencode/repeats/ChromHMM/
#                        SPIN/GC --, and (stages 12-13, only if histone_marks.txt/
#                        speckle_names.txt are present) a metagene histone-mark/
#                        speckle-track correlation report; independent of stage 07, runs
#                        alongside it)
#
# Every stage writes into the CURRENT directory (matching each stage script's own
# convention of writing next to its inputs) -- run this from wherever you want the
# pipeline's output to land, and make sure config paths point at your real data.
#
# LSF parallelization: every stage that loops over independent items (samples / batches /
# constructs / COMBOS rows -- each reads its own inputs and writes its own outputs, no
# shared state) submits one bsub job per item instead of looping sequentially, then blocks
# until that whole batch finishes before the next (dependent) stage starts -- so this
# script itself stays a single long-lived foreground process (submit it via bsub yourself,
# or nohup/tmux it, if you want to disconnect), but the actual per-item work runs in
# parallel across the cluster. Stages with no independent-item loop (e.g. stage 11, which
# processes every construct in one Rscript call) still get submitted as a single LSF job
# rather than run directly on the login node. Set USE_LSF=0 in the environment to fall back
# to plain sequential local execution instead (e.g. for a small test run without LSF
# access).
# =============================================================================
set -euo pipefail

CONFIG=${1:?"usage: run_pipeline.sh <config.txt> [target] [threads]"}
TARGET=${2:-combined_counts}
THREADS=${3:-8}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"
WORK=$(pwd)

case "$TARGET" in
  cutadapt|mapping|genuine_barcodes|combined_counts|annotated) ;;
  *) echo "unknown target '$TARGET' (expected: cutadapt | mapping | genuine_barcodes | combined_counts | annotated)" >&2; exit 1;;
esac

echo "=== TRIP unified pipeline ==="
echo "config: $CONFIG"
echo "target: $TARGET"
echo "workdir: $WORK"
echo

# -----------------------------------------------------------------------------
# Parse config into section files under a scratch dir, and source ### GLOBAL as env vars.
# -----------------------------------------------------------------------------
PARSED=$(mktemp -d)
trap 'rm -rf "$PARSED"' EXIT

for section in GLOBAL FASTQ_MAPS FASTQ_NORM FASTQ_SAMPLES COMBOS SPLICED_UNSPLICED; do
  awk -v s="### $section" '
    $0==s {grab=1; next}
    /^### /{grab=0}
    grab && NF && $0 !~ /^#/ {print}
  ' "$CONFIG" > "$PARSED/$section.tsv"
done

# GLOBAL: key<TAB>value -> environment variables
while IFS=$'\t' read -r key value; do
  [ -z "$key" ] && continue
  export "$key"="$value"
done < "$PARSED/GLOBAL.tsv"

: "${BARCODE_LEN:?BARCODE_LEN missing from ### GLOBAL}"
: "${BOWTIE2_INDEX:?BOWTIE2_INDEX missing from ### GLOBAL}"
: "${CONT_FACTOR:=10}"
: "${R1_PRIMER:?R1_PRIMER missing from ### GLOBAL}"
: "${R1_SUFFIX:?R1_SUFFIX missing from ### GLOBAL}"
: "${R2_PRIMER:?R2_PRIMER missing from ### GLOBAL}"
TTAA_BED=${TTAA_BED:-}
DPNII_BED=${DPNII_BED:-}
# Stages 09-11 (annotated target only) -- each left blank falls back to that stage
# script's own hardcoded default (same convention as DPNII_BED above).
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
SPECKLE_DIR=${SPECKLE_DIR:-}
HISTONE_BW_DIR=${HISTONE_BW_DIR:-}
HG38_TO_HG19_CHAIN=${HG38_TO_HG19_CHAIN:-}

echo "--- GLOBAL settings ---"
echo "BARCODE_LEN=$BARCODE_LEN  CONT_FACTOR=$CONT_FACTOR  TTAA_BED=${TTAA_BED:-<off>}  DPNII_BED=${DPNII_BED:-<script default>}"
echo "BOWTIE2_INDEX=$BOWTIE2_INDEX"
if [ "$TARGET" = "annotated" ]; then
  echo "MANE_BED=${MANE_BED:-<script default>}  MANE_INTRONS_BED=${MANE_INTRONS_BED:-<script default>}  GENCODE_NONCODING_BED=${GENCODE_NONCODING_BED:-<script default>}"
  echo "REPEATS_BED=${REPEATS_BED:-<script default>}  CHROM_SIZES=${CHROM_SIZES:-<script default>}  GENOME_FASTA=${GENOME_FASTA:-<script default>}"
  echo "TRANS2GENE=${TRANS2GENE:-<script default>}  MANE_EXPRESSION=${MANE_EXPRESSION:-<script default>}  CHROMHMM_BED=${CHROMHMM_BED:-<script default>}  SPIN_BED=${SPIN_BED:-<script default>}"
  echo "SPECKLE_DIR=${SPECKLE_DIR:-<script default>}  HISTONE_BW_DIR=${HISTONE_BW_DIR:-<script default>}  HG38_TO_HG19_CHAIN=${HG38_TO_HG19_CHAIN:-<script default>}"
fi
echo

export R1_PRIMER R1_SUFFIX R2_PRIMER
export BARCODE_LEN

# -----------------------------------------------------------------------------
# LSF parallel-submission helpers -- every stage below that loops over independent items
# (samples / batches / constructs / COMBOS rows -- each reads its own inputs and writes its
# own outputs, no cross-item state) submits one bsub job per item instead of running them in
# a sequential loop, then blocks until the whole batch finishes before the next
# (dependent) stage starts. USE_LSF=0 falls back to plain sequential execution (e.g. for a
# small local/test run without LSF access) -- default is on, matching "parallelize
# wherever it helps".
# -----------------------------------------------------------------------------
USE_LSF=${USE_LSF:-1}
# LSF's own configured default queue is "short"; this cluster's work belongs on
# gsla-cpu, so pass it explicitly rather than relying on the site default.
LSF_QUEUE=${LSF_QUEUE:-gsla-cpu}
declare -a JOBIDS=()
declare -A JOBLOG=()

# Each bsub job gets a fresh, non-interactive shell that does NOT inherit whatever lmod
# modules happen to be loaded in the shell that launched this script, so every command
# below loads what it needs explicitly (module names below match what's actually available
# via `module avail` on this cluster; bump these if the cluster's module set changes).
CUTADAPT_MOD="cutadapt/4.2-GCCcore-11.3.0"
SAMTOOLS_MOD="SAMtools/1.21-GCC-13.2.0"
# 02_run_bowtie2.sh also pipes through samtools view/sort/index and bedtools bamtobed (see
# that script's own header -- it's a reconstruction of an ad hoc alignment step that never
# ran through this orchestrator before), so its module list needs to cover all three, not
# just Bowtie2 itself -- caught 2026-08-07 when stage 02 failed for real for the first time
# (SON_run4) with "samtools: command not found".
BOWTIE2_MOD="Bowtie2/2.5.1-GCC-12.3.0 $SAMTOOLS_MOD BEDTools/2.31.1-GCC-13.3.0"
R_MOD="R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1"
BEDTOOLS_MOD="BEDTools/2.31.1-GCC-13.3.0"
KENTUTILS_MOD="kentUtils/v377"

# Extra bsub select clause for the stages that load a Bioconductor package with compiled
# code (rtracklayer, ComplexHeatmap). model==AMD_EPYC below is NOT a fine enough filter: it
# matches both amd_epyc_9k_* (EPYC 9004, Zen 4, has AVX-512) and amd_epyc_7k_* (cn650-cn653,
# EPYC 7702, Zen 2, does NOT). R-bundle-Bioconductor/3.19's rtracklayer .so carries AVX-512
# instructions, so on a 7k host dyn.load() dies instantly with "*** caught illegal operation
# *** cause 'illegal operand'" and LSF exit code 132 (SIGILL) -- observed 2026-08-19, stage 11
# on cn651. ncpus>=160 is the discriminator: every 9k CPU host has 168 (84c) or 256 (128c)
# ncpus while every 7k host has 128 or fewer, and unlike `bsub -m amd_epyc_9k_*_cpu_hosts`
# (rejected with "Host or host group is not used by the queue") it works from gsla-cpu, whose
# only host group is gsla_cpu_hosts. The pure-R and bedtools stages run fine on either
# generation, so only the Bioconductor stages pay the narrower selection.
BIOC_SELECT=${BIOC_SELECT:-" && ncpus>=160"}

# submitOrRun <logfile> <mem_mb> <n_cores> <modules|-> <command...>
# <modules> is a space-separated list of lmod module names to load before running <command>
# (pass "-" for none). Applied in both the bsub and local-eval branches so USE_LSF=0 behaves
# the same way (harmless if the module is already loaded interactively).
# Prefix a call with SUBMIT_SELECT="<extra select terms>" (see BIOC_SELECT above) to narrow
# that one job's host selection; unset, the job goes to any host matching -R below.
submitOrRun() {
  local logfile=$1 mem=$2 ncores=$3 modules=$4
  shift 4
  local preamble=""
  if [ "$modules" != "-" ] && [ -n "$modules" ]; then
    preamble="source /etc/profile.d/lmod.sh; module load $modules; "
  fi
  if [ "$USE_LSF" = "1" ]; then
    local out jobid
    # bsub -o APPENDS to an existing file (only -oo overwrites), so a stale log from an
    # earlier run would still carry its old "Successfully completed."/"Exited with exit
    # code" epilogue -- which is exactly what waitForJobs reads to decide this job's fate.
    # Remove it first so every marker found belongs to THIS submission.
    rm -f "$logfile"
    # Compiled-binary modules (bedtools/bowtie2/cutadapt/kentUtils) live under an
    # AMD-specific easybuild tree and aren't visible on non-AMD_EPYC hosts, so every job is
    # pinned to that host class -- cheap and correct even for the pure-R stages.
    out=$(bsub -q "$LSF_QUEUE" -R "select[type==X86_64 && model==AMD_EPYC${SUBMIT_SELECT:-}] rusage[mem=${mem}]" -n "$ncores" -o "$logfile" "${preamble}$*" 2>&1)
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
# exit code N" / "TERM_*"). Aborts the whole pipeline (set -e via exit 1) if any job failed,
# since a silently-missing output would just surface as a confusing failure several stages
# later.
#
# Deliberately does NOT decide anything from `bjobs -o stat <id>`, which is what this
# function used to do and what aborted a real IntronTRIP run on 2026-08-04 after all 32 of
# its stage-04 jobs had in fact succeeded: `bjobs <id>` yields an EMPTY status both once LSF
# has purged the finished job's record and whenever mbatchd is transiently unreachable
# ("mbatchd is down or the ports are inaccessable", seen on this cluster the same morning)
# -- and empty was being mapped to UNKNOWN, i.e. failure. Neither condition says anything
# about whether the job worked. The job's own log is authoritative and never expires;
# `bjobs -w` (active listing, RUN/PEND only) is used purely as the "is it still going?"
# signal, never as the verdict.
waitForJobs() {
  [ "$USE_LSF" = "1" ] || { JOBIDS=(); return 0; }
  [ ${#JOBIDS[@]} -eq 0 ] && return 0
  local failed=0 jid log active gone

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
    echo "One or more jobs in this stage failed -- aborting pipeline." >&2
    exit 1
  fi
}

# -----------------------------------------------------------------------------
# Stage 01/02 helper: link a sample's raw fastqs into place under path_prefix's basename,
# since 01_run_cutadapt.sh expects "<name>_R1_001.fastq.gz" etc in the CURRENT directory.
# -----------------------------------------------------------------------------
stageFastq() {
  local name=$1 pathPrefix=$2 r1suf=$3 r2suf=$4
  local r1src="${pathPrefix}${r1suf}" r2src="${pathPrefix}${r2suf}"
  [ -f "$r1src" ] || { echo "missing R1 for $name: $r1src" >&2; exit 1; }
  [ -f "$r2src" ] || { echo "missing R2 for $name: $r2src" >&2; exit 1; }
  ln -sf "$r1src" "${name}_R1_001.fastq.gz"
  ln -sf "$r2src" "${name}_R2_001.fastq.gz"
}

# =============================================================================
# Stage 01 (+ barcode/count extraction): ### FASTQ_MAPS in map mode, ### FASTQ_NORM + ### FASTQ_SAMPLES
# in count mode.
# =============================================================================
echo "=== Stage 01: cutadapt (one job per sample) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name (mode=map) --"
  stageFastq "$name" "$pathPrefix" "$r1suf" "$r2suf"
  submitOrRun "${name}_01_cutadapt.log" 8000 "$THREADS" "$CUTADAPT_MOD" bash "$SCRIPTS/01_run_cutadapt.sh" map "$name" "$THREADS"
done < "$PARSED/FASTQ_MAPS.tsv"

for section in FASTQ_NORM FASTQ_SAMPLES; do
  while IFS=$'\t' read -r name pathPrefix r1suf r2suf; do
    [ -z "$name" ] && continue
    echo "-- $name (mode=count) --"
    stageFastq "$name" "$pathPrefix" "$r1suf" "$r2suf"
    # Optional per-row 5th/6th columns: R1/R2 adapter override (blank = stage 01's
    # default). Looked up via awk rather than folded into the `read` above -- bash's
    # `read` collapses a genuinely empty field between two tabs regardless of IFS
    # (tab is always "IFS whitespace"), which would silently shift r2Adapter's value
    # into r1Adapter whenever r1Adapter is blank and r2Adapter isn't -- exactly the
    # common case here (e.g. every IntronTRIP_combined/SON_combined Spliced/Unspliced
    # row). awk's field splitting doesn't have this collapsing behavior.
    r1Adapter=$(awk -F'\t' -v n="$name" '$1==n{print $5}' "$PARSED/$section.tsv")
    r2Adapter=$(awk -F'\t' -v n="$name" '$1==n{print $6}' "$PARSED/$section.tsv")
    submitOrRun "${name}_01_cutadapt.log" 8000 "$THREADS" "$CUTADAPT_MOD" \
      env R1_ADAPTER="$r1Adapter" R2_ADAPTER="$r2Adapter" bash "$SCRIPTS/01_run_cutadapt.sh" count "$name" "$THREADS"
  done < "$PARSED/$section.tsv"
done
waitForJobs
echo

if [ "$TARGET" = "cutadapt" ]; then echo "target=cutadapt reached, stopping."; exit 0; fi

echo "=== Stage 02: bowtie2 (### FASTQ_MAPS samples only, one job per sample) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name --"
  submitOrRun "${name}_02_bowtie2.log" 16000 "$THREADS" "$BOWTIE2_MOD" bash "$SCRIPTS/02_run_bowtie2.sh" "$name" "$THREADS" "$BOWTIE2_INDEX"
done < "$PARSED/FASTQ_MAPS.tsv"
waitForJobs
echo

if [ "$TARGET" = "mapping" ]; then echo "target=mapping reached, stopping."; exit 0; fi

# =============================================================================
# Stage 03/04: per ### FASTQ_MAPS sample
# =============================================================================
echo "=== Stage 03: analyzeTRIPMapping (one job per sample) ==="
# 03_analyzeTRIPMapping.R reads its two inputs with base read.table() rather than
# data.table::fread(), which is far more memory-hungry per row -- 8000 (the original
# default, also never stress-tested for real before) reliably OOM'd (TERM_MEMLIMIT) on all
# 4 SON_run4 Map libraries despite each decompressing to well under 1GB (~6-7M lines each).
# Bumped generously rather than tuned tight, matching this project's usual OOM-recovery
# margin for base-R barcode-table jobs.
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name --"
  submitOrRun "${name}_03_analyzeTRIPMapping.log" 32000 1 "$R_MOD" \
    Rscript "$SCRIPTS/03_analyzeTRIPMapping.R" "$WORK/" "$WORK/02_mapping/" "$name" "$BARCODE_LEN"
done < "$PARSED/FASTQ_MAPS.tsv"
waitForJobs
echo

echo "=== Stage 04: filterMappings (+ TTAA if configured; one job per sample) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name --"
  if [ -n "$TTAA_BED" ]; then
    submitOrRun "${name}_04_filterMappings.log" 4000 1 "$R_MOD" \
      Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$WORK/" "$name" ttaaFilter "$TTAA_BED"
  else
    submitOrRun "${name}_04_filterMappings.log" 4000 1 "$R_MOD" \
      Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$WORK/" "$name"
  fi
done < "$PARSED/FASTQ_MAPS.tsv"
waitForJobs
echo

# =============================================================================
# Stage 05: filterContaminants -- one run PER contamination batch (### FASTQ_MAPS' 6th column).
# Index-hopping contamination is a sequencing-run/flowcell artefact: samples that were
# never multiplexed together can't contaminate each other, so real historical fileLists
# (e.g. run3.txt/run5.txt, IntronTRIP_run2.txt/_run3.txt/_run4.txt) are scoped one per
# physical run, not one global fileList across every map sample in a combined dataset.
# Samples with no batch value fall into a single "default" batch. The contamination group
# for each sample is its mapName (### FASTQ_MAPS' 5th column) -- the two concepts are the same
# thing here, so there's no separate "group" column.
# =============================================================================
echo "=== Stage 05: filterContaminants (one job per batch) ==="
BATCHES=$(awk -F'\t' '{print ($6==""?"default":$6)}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
if [ -z "$BATCHES" ]; then
  echo "no map samples, skipping"
else
  while read -r batch; do
    [ -z "$batch" ] && continue
    FILELIST="$WORK/.pipeline_fileList.${batch}.txt"
    awk -F'\t' -v b="$batch" '($6==""?"default":$6)==b {print $1"\t"$5}' "$PARSED/FASTQ_MAPS.tsv" > "$FILELIST"
    echo "-- batch: $batch ($(wc -l < "$FILELIST") sample(s)) --"
    submitOrRun "stage05_${batch}.log" 8000 1 "$R_MOD" \
      Rscript "$SCRIPTS/05_filterContaminants_uniform.R" "$WORK/" "$(basename "$FILELIST")" "$CONT_FACTOR"
  done <<< "$BATCHES"
fi
waitForJobs
echo

# =============================================================================
# Stage 06: get_genuine_barcodes -- one run per distinct mapName in ### FASTQ_MAPS
# (one-to-many: every sample sharing a mapName is combined into that construct's map)
# =============================================================================
echo "=== Stage 06: get_genuine_barcodes (one job per construct) ==="
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
while read -r mapName; do
  [ -z "$mapName" ] && continue
  manifest="$WORK/${mapName}.txt"
  awk -F'\t' -v m="$mapName" '$5==m{print $1".retained_barcodes.bed"}' "$PARSED/FASTQ_MAPS.tsv" > "$manifest"
  echo "-- $mapName ($(wc -l < "$manifest") sample(s)) --"
  submitOrRun "stage06_${mapName}.log" 8000 1 "$R_MOD" \
    Rscript "$SCRIPTS/06_get_genuine_barcodes.R" "$WORK/" "$(basename "$manifest")" "$BARCODE_LEN" "$DPNII_BED"
done <<< "$MAPNAMES"
waitForJobs
echo

if [ "$TARGET" = "genuine_barcodes" ]; then echo "target=genuine_barcodes reached, stopping."; exit 0; fi

# =============================================================================
# Stage 07: normalizeTRIP -- one run per ### COMBOS row
# =============================================================================
# One job per COMBOS row -- each row reads its own map/norm/exp inputs and writes its own
# <expSample>.combined_counts.txt, no state shared between rows.
echo "=== Stage 07: normalizeTRIP (combined_counts.txt, one job per COMBOS row) ==="
while IFS=$'\t' read -r mapName normSample expSample effLength revComp oneMismatch; do
  [ -z "$mapName" ] && continue
  echo "-- map=$mapName norm=$normSample exp=$expSample --"
  args=("$WORK/" "$WORK/" "${mapName}.txt" "$normSample" "$expSample" "$effLength" "$revComp")
  if [ "${oneMismatch:-}" = "yes" ]; then args+=("oneMismatch"); fi
  submitOrRun "stage07_${expSample}.log" 8000 1 "$R_MOD" Rscript "$SCRIPTS/07_normalizeTRIP_MM.R" "${args[@]}"
done < "$PARSED/COMBOS.tsv"
waitForJobs

echo
echo "=== combined_counts.txt written for each ### COMBOS row ==="

# =============================================================================
# Stage 08: combineSplicedUnsplicedPairs -- optional, only if ### SPLICED_UNSPLICED has
# rows. Depends on stage 07's combined_counts.txt (the spliced/unspliced/total triplets
# it names), so it runs for both "combined_counts" and "annotated" targets, right after
# stage 07 -- unlike stages 09-13, which only need stage 06 and are annotated-only.
# =============================================================================
if [ -s "$PARSED/SPLICED_UNSPLICED.tsv" ]; then
  echo "=== Stage 08: combineSplicedUnsplicedPairs (spliced/unspliced/total ratio files) ==="
  cp "$PARSED/SPLICED_UNSPLICED.tsv" "$WORK/spliced_unspliced.txt"
  # Single job (this script loops over every ### SPLICED_UNSPLICED triplet internally), but
  # still submitted rather than run inline -- it used to be the one stage that called Rscript
  # directly, which contradicted this script's own "no stage runs on the login node" contract
  # above and meant it silently skipped the module load every other stage gets.
  submitOrRun "stage08_combineSplicedUnspliced.log" 8000 1 "$R_MOD" \
    Rscript "$SCRIPTS/08_combineSplicedUnsplicedPairs_uniform.R" "$WORK/"
  waitForJobs
  echo
fi

if [ "$TARGET" != "annotated" ]; then echo "target=combined_counts reached, stopping."; exit 0; fi

# =============================================================================
# Stages 09-11 (target=annotated): per-construct reference annotation, independent of
# stage 07 -- only needs stage 06's genuine-barcode maps, not combined_counts.txt.
# Reuses $MAPNAMES from stage 06. Builds one name-list file (one "<mapName>.txt" line
# per construct, matching the real production names.complete_maps.txt format) that
# stages 09 and 11 each take as a single whole-list argument; stage 10 runs once per
# construct instead, matching its original per-construct bsub usage.
# =============================================================================
NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

echo
# Stage 09 loops internally (8 sequential intersects x every construct) inside its own
# script rather than per-construct in this orchestrator, so it isn't split into one bsub
# job per item the way the other stages are -- would need 09_make_tracks_genuine_barcodes.sh
# itself restructured to take a single construct name to parallelize further. Still routed
# through submitOrRun so it runs as one LSF job rather than blocking on the login node.
echo "=== Stage 09: make_tracks (MANE/gencode/repeats intersects) ==="
submitOrRun "stage09_make_tracks.log" 4000 1 "$BEDTOOLS_MOD" bash "$SCRIPTS/09_make_tracks_genuine_barcodes.sh" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$REPEATS_BED"
waitForJobs
echo

echo "=== Stage 10: gc_content (one job per construct) ==="
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "-- $mapName --"
  submitOrRun "stage10_${mapName}.log" 4000 1 "$BEDTOOLS_MOD" bash "$SCRIPTS/10_gc_content_single_genuine_barcodes.sh" "${mapName}.txt" "$CHROM_SIZES" "$GENOME_FASTA"
done <<< "$MAPNAMES"
waitForJobs
echo

echo "=== Stage 11: table_annotation (annotated_tables.txt + annotated_table_used.txt) ==="
NORMMANIFEST="$WORK/.pipeline_norm_manifest.txt"
awk -F'\t' '{print $1"\t"$2}' "$PARSED/COMBOS.tsv" > "$NORMMANIFEST"
SUBMIT_SELECT="$BIOC_SELECT" submitOrRun "stage11_table_annotation.log" 16000 1 "$R_MOD" Rscript "$SCRIPTS/11_table_annotation_genuine_barcodes.R" "$WORK/" "$(basename "$NAMELIST")" \
  "$MANE_BED" "$MANE_INTRONS_BED" "$GENCODE_NONCODING_BED" "$TRANS2GENE" "$MANE_EXPRESSION" "$CHROMHMM_BED" "$SPIN_BED" "$(basename "$NORMMANIFEST")"
waitForJobs
echo

# =============================================================================
# comparePairwiseSharedPositions / annotateTRIP_positions / annotateTRIPGeneBodies /
# gc_content_report / compareReplicatesAndMultiIntegrations (previously numbered as pipeline stages, now unnumbered)
# are visualization/QC scripts, not pipeline stages -- run them by hand against stage
# 07/10's output (see README's "Visualization scripts" section), same as
# plotSplicingEfficiencies_uniform.R and computeIntegrationStats.R. Not orchestrated here.
# =============================================================================

# =============================================================================
# Stages 12-13: metagene bigWig overlap (histone marks + speckle/TSA-seq tracks) and its
# correlation-vs-expression report. Optional: skipped entirely if histone_marks.txt and
# speckle_names.txt (real, curated per-dataset files -- see e.g.
# single_combined/histone_marks.txt -- not generated by this pipeline) aren't both present
# in the working directory. Builds its own <mapName><TAB><expName> $EXPLIST from ### COMBOS'
# mapName+expSample columns (1st+3rd) -- stage 13's expFile is that same mapType<TAB>exp
# format.
# =============================================================================
EXPLIST="$WORK/.pipeline_exp_list.txt"
awk -F'\t' '{print $1"\t"$3}' "$PARSED/COMBOS.tsv" > "$EXPLIST"
ALLMAPNAMES=$(echo "$MAPNAMES" | paste -sd,)
if [ -f "$WORK/histone_marks.txt" ] && [ -f "$WORK/speckle_names.txt" ]; then
  echo "=== Stage 12: meta_bigwig_overlap (bedtools slop + bigWigAverageOverBed, one job per construct) ==="
  # 12_meta_bigwig_overlap.sh's own header flags this: it's a straight sequential port of
  # the original csh scripts, looping over every (construct x histone/speckle track) pair
  # as one bigWigAverageOverBed call after another inside a single job -- e.g. 4
  # constructs x 59 tracks (single_combined) is 236 sequential external-tool calls in one
  # job. The script itself already accepts a MAPTYPES_FILE with just one mapName, so no
  # change to its own logic is needed -- fan out one job per construct here instead of a
  # single job over the whole list, matching stage 10 (gc_content)'s existing
  # one-job-per-construct convention just above.
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
  SUBMIT_SELECT="$BIOC_SELECT" submitOrRun "stage13_metagene.log" 8000 1 "$R_MOD" Rscript "$SCRIPTS/13_plotMetageneCorrelations_uniform.R" "$WORK/" pearson "$ALLMAPNAMES" \
    "$(basename "$EXPLIST")" metagene_correlations
  waitForJobs
  echo
else
  echo "=== Stages 12-13 skipped: histone_marks.txt and/or speckle_names.txt not found in $WORK ==="
  echo
fi

echo "=== done: <mapName>.txt.annotated_tables.txt and .annotated_table_used.txt written for each construct, and (if histone_marks.txt/speckle_names.txt present) metagene_correlations.{meta_corrs,summary_meta_corrs}.pdf written ==="
