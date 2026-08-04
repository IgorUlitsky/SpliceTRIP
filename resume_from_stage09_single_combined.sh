#!/bin/bash
# One-off resume script: re-run single_combined stages 09-11 (the core data-producing
# stages: MANE/gencode/repeat intersects, GC content, annotated_tables.txt), skipping 01-07.
# Needed because stage 06's retained_barcodes_genuine.bed was regenerated (MANE_intronic fix
# + position-proximity dedup + TTAA filter + contamination fix) but stages 09-11 were never
# re-run after that -- confirmed via mtimes: stage09's MANE_all.txt etc are from 2026-07-24,
# stage10's _with_gc.bed from 2026-07-27, both predating stage06's 2026-08-03 regeneration.
# The visualization/QC scripts (unnumbered) are (PDFs, mostly) built on top of annotated_tables.txt --
# not rerun here, can be refreshed separately if the plots themselves are needed.
#
# Mirrors run_pipeline.sh's own stage 09-11 block (including its LSF parallel-submission
# pattern) rather than reimplementing it differently, minus stages 01-07's setup.
set -euo pipefail

CONFIG=${1:?"usage: resume_from_stage09_single_combined.sh <config.txt> <single_combined_dir>"}
WORK=${2:?"usage: resume_from_stage09_single_combined.sh <config.txt> <single_combined_dir>"}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"

echo "=== Resume from stage 09 ==="
echo "config: $CONFIG"
echo "single_combined dir: $WORK"
echo

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

USE_LSF=${USE_LSF:-1}
declare -a JOBIDS=()

submitOrRun() {
  local logfile=$1 mem=$2 ncores=$3
  shift 3
  if [ "$USE_LSF" = "1" ]; then
    local out jobid
    # BEDTools/R modules live under an AMD-specific easybuild tree, invisible on non-AMD_EPYC
    # hosts -- pin every job there (cheap, and correct for both the bedtools and R stages).
    out=$(bsub -R "select[type==X86_64 && model==AMD_EPYC] rusage[mem=${mem}]" -n "$ncores" -o "$logfile" "$@" 2>&1)
    jobid=$(echo "$out" | grep -oE '<[0-9]+>' | head -1 | tr -d '<>')
    if [ -z "$jobid" ]; then
      echo "ERROR: bsub submission failed:" >&2
      echo "$out" >&2
      exit 1
    fi
    JOBIDS+=("$jobid")
  else
    eval "$@" > "$logfile" 2>&1
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
    echo "One or more jobs in this stage failed -- aborting." >&2
    exit 1
  fi
}

cd "$WORK"
MAPNAMES=$(awk -F'\t' '{print $5}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)

NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

# bsub jobs get a fresh shell that doesn't inherit interactively-loaded lmod modules, so each
# command loads what it needs explicitly (same pattern proven working in
# resume_stage07_parallel_single_combined.sh's direct bsub calls).
LMOD_INIT="source /etc/profile.d/lmod.sh"
BEDTOOLS_MOD="BEDTools/2.31.1-GCC-13.3.0"
R_MOD="R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1"

echo "=== Stage 09: make_tracks (MANE/gencode/repeats intersects) ==="
submitOrRun "stage09_make_tracks.log" 4000 1 \
  "$LMOD_INIT; module load $BEDTOOLS_MOD; bash $SCRIPTS/09_make_tracks_genuine_barcodes.sh $(basename "$NAMELIST") $MANE_BED $MANE_INTRONS_BED $GENCODE_NONCODING_BED $REPEATS_BED"
waitForJobs
echo

echo "=== Stage 10: gc_content (one job per construct) ==="
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "-- $mapName --"
  submitOrRun "stage10_${mapName}.log" 4000 1 \
    "$LMOD_INIT; module load $BEDTOOLS_MOD; bash $SCRIPTS/10_gc_content_single_genuine_barcodes.sh ${mapName}.txt $CHROM_SIZES $GENOME_FASTA"
done <<< "$MAPNAMES"
waitForJobs
echo

echo "=== Stage 11: table_annotation (annotated_tables.txt + annotated_table_used.txt) ==="
NORMMANIFEST="$WORK/.pipeline_norm_manifest.txt"
awk -F'\t' '{print $1"\t"$2}' "$PARSED/COMBOS.tsv" > "$NORMMANIFEST"
submitOrRun "stage11_table_annotation.log" 16000 1 \
  "$LMOD_INIT; module load $R_MOD; Rscript $SCRIPTS/11_table_annotation_genuine_barcodes.R $WORK/ $(basename "$NAMELIST") $MANE_BED $MANE_INTRONS_BED $GENCODE_NONCODING_BED $TRANS2GENE $MANE_EXPRESSION $CHROMHMM_BED $SPIN_BED $(basename "$NORMMANIFEST")"
waitForJobs
echo

echo "=== done (stages 09-11 only; the unnumbered visualization/QC scripts are not rerun here) ==="
