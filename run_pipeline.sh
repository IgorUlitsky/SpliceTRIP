#!/bin/bash
# =============================================================================
# TRIP unified pipeline — orchestrator
# =============================================================================
# Reads a config text file (see config/example_config.txt) and runs stages 01-11
# (cutadapt -> bowtie2 -> analyzeTRIPMapping -> filterMappings[+TTAA] ->
# filterContaminants -> get_genuine_barcodes -> normalizeTRIP[+oneMismatch] ->
# make_tracks -> gc_content -> table_annotation -> comparePairwiseSharedPositions) in order.
#
# Usage
#   bash run_pipeline.sh <config.txt> [target] [threads]
#
#   target (default: combined_counts) — how far to run:
#     cutadapt          stage 01 only
#     mapping           01-02 (Map-mode samples only)
#     genuine_barcodes  01-06
#     combined_counts   01-07 (everything needed for combined_counts.txt) ★ pipeline's main deliverable
#     annotated         01-11 (adds per-construct annotated_tables.txt -- MANE/gencode/
#                        repeats/ChromHMM/SPIN/GC -- and a pairwise shared-position
#                        comparison across every ### FASTQ_MAPS construct; independent of
#                        stage 07, runs alongside it)
#
# Every stage writes into the CURRENT directory (matching each stage script's own
# convention of writing next to its inputs) -- run this from wherever you want the
# pipeline's output to land, and make sure config paths point at your real data.
#
# This runs everything locally/sequentially. In production, stages 01-02 (cutadapt,
# bowtie2) are heavy (tens of minutes, several GB) and should be submitted as LSF jobs
# instead of run directly on a login node -- this script does not submit jobs itself;
# wrap the relevant section in bsub yourself, or ask for that to be added once this is
# validated locally.
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

for section in GLOBAL FASTQ_MAPS FASTQ_NORM FASTQ_SAMPLES COMBOS; do
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
# Stages 08-10 (annotated target only) -- each left blank falls back to that stage
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

echo "--- GLOBAL settings ---"
echo "BARCODE_LEN=$BARCODE_LEN  CONT_FACTOR=$CONT_FACTOR  TTAA_BED=${TTAA_BED:-<off>}  DPNII_BED=${DPNII_BED:-<script default>}"
echo "BOWTIE2_INDEX=$BOWTIE2_INDEX"
if [ "$TARGET" = "annotated" ]; then
  echo "MANE_BED=${MANE_BED:-<script default>}  MANE_INTRONS_BED=${MANE_INTRONS_BED:-<script default>}  GENCODE_NONCODING_BED=${GENCODE_NONCODING_BED:-<script default>}"
  echo "REPEATS_BED=${REPEATS_BED:-<script default>}  CHROM_SIZES=${CHROM_SIZES:-<script default>}  GENOME_FASTA=${GENOME_FASTA:-<script default>}"
  echo "TRANS2GENE=${TRANS2GENE:-<script default>}  MANE_EXPRESSION=${MANE_EXPRESSION:-<script default>}  CHROMHMM_BED=${CHROMHMM_BED:-<script default>}  SPIN_BED=${SPIN_BED:-<script default>}"
fi
echo

export R1_PRIMER R1_SUFFIX R2_PRIMER
export BARCODE_LEN

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
echo "=== Stage 01: cutadapt ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name (mode=map) --"
  stageFastq "$name" "$pathPrefix" "$r1suf" "$r2suf"
  bash "$SCRIPTS/01_run_cutadapt.sh" map "$name" "$THREADS"
done < "$PARSED/FASTQ_MAPS.tsv"

for section in FASTQ_NORM FASTQ_SAMPLES; do
  while IFS=$'\t' read -r name pathPrefix r1suf r2suf; do
    [ -z "$name" ] && continue
    echo "-- $name (mode=count) --"
    stageFastq "$name" "$pathPrefix" "$r1suf" "$r2suf"
    bash "$SCRIPTS/01_run_cutadapt.sh" count "$name" "$THREADS"
  done < "$PARSED/$section.tsv"
done
echo

if [ "$TARGET" = "cutadapt" ]; then echo "target=cutadapt reached, stopping."; exit 0; fi

echo "=== Stage 02: bowtie2 (### FASTQ_MAPS samples only) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name --"
  bash "$SCRIPTS/02_run_bowtie2.sh" "$name" "$THREADS" "$BOWTIE2_INDEX"
done < "$PARSED/FASTQ_MAPS.tsv"
echo

if [ "$TARGET" = "mapping" ]; then echo "target=mapping reached, stopping."; exit 0; fi

# =============================================================================
# Stage 03/04: per ### FASTQ_MAPS sample
# =============================================================================
echo "=== Stage 03: analyzeTRIPMapping ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name --"
  Rscript "$SCRIPTS/03_analyzeTRIPMapping.R" "$WORK/" "$WORK/02_mapping/" "$name" "$BARCODE_LEN"
done < "$PARSED/FASTQ_MAPS.tsv"
echo

echo "=== Stage 04: filterMappings (+ TTAA if configured) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mapName batch; do
  [ -z "$name" ] && continue
  echo "-- $name --"
  if [ -n "$TTAA_BED" ]; then
    Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$WORK/" "$name" ttaaFilter "$TTAA_BED"
  else
    Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$WORK/" "$name"
  fi
done < "$PARSED/FASTQ_MAPS.tsv"
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
echo "=== Stage 05: filterContaminants ==="
BATCHES=$(awk -F'\t' '{print ($6==""?"default":$6)}' "$PARSED/FASTQ_MAPS.tsv" | sort -u)
if [ -z "$BATCHES" ]; then
  echo "no map samples, skipping"
else
  while read -r batch; do
    [ -z "$batch" ] && continue
    FILELIST="$WORK/.pipeline_fileList.${batch}.txt"
    awk -F'\t' -v b="$batch" '($6==""?"default":$6)==b {print $1"\t"$5}' "$PARSED/FASTQ_MAPS.tsv" > "$FILELIST"
    echo "-- batch: $batch ($(wc -l < "$FILELIST") sample(s)) --"
    Rscript "$SCRIPTS/05_filterContaminants_uniform.R" "$WORK/" "$(basename "$FILELIST")" "$CONT_FACTOR"
  done <<< "$BATCHES"
fi
echo

# =============================================================================
# Stage 06: get_genuine_barcodes -- one run per distinct mapName in ### FASTQ_MAPS
# (one-to-many: every sample sharing a mapName is combined into that construct's map)
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

if [ "$TARGET" = "genuine_barcodes" ]; then echo "target=genuine_barcodes reached, stopping."; exit 0; fi

# =============================================================================
# Stage 07: normalizeTRIP -- one run per ### COMBOS row
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
echo "=== combined_counts.txt written for each ### COMBOS row ==="

if [ "$TARGET" != "annotated" ]; then echo "target=combined_counts reached, stopping."; exit 0; fi

# =============================================================================
# Stages 08-10 (target=annotated): per-construct reference annotation, independent of
# stage 07 -- only needs stage 06's genuine-barcode maps, not combined_counts.txt.
# Reuses $MAPNAMES from stage 06. Builds one name-list file (one "<mapName>.txt" line
# per construct, matching the real production names.complete_maps.txt format) that
# stages 08 and 10 each take as a single whole-list argument; stage 09 runs once per
# construct instead, matching its original per-construct bsub usage.
# =============================================================================
NAMELIST="$WORK/.pipeline_names_complete_maps.txt"
> "$NAMELIST"
while read -r mapName; do
  [ -z "$mapName" ] && continue
  echo "${mapName}.txt" >> "$NAMELIST"
done <<< "$MAPNAMES"

echo
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

# =============================================================================
# Stage 11: comparePairwiseSharedPositions -- every pair of ### FASTQ_MAPS mapNames (not
# every ### COMBOS row), so this needs its own name-list and its own <mapName><TAB>
# <expName> exp-list, built from ### COMBOS' mapName+expSample columns (1st+3rd). Real
# production expListFiles are often hand-curated and NOT 1:1 with ### COMBOS (e.g. a
# separate splicing-ratio pipeline's outputs pooled in) -- see examples/exp_lists/ for
# real ones to use instead of this auto-generated one, by passing your own expListFile as
# a 3rd arg to this script directly rather than going through run_pipeline.sh.
# =============================================================================
echo "=== Stage 11: comparePairwiseSharedPositions ==="
mkdir -p "$WORK/plots"
EXPLIST="$WORK/.pipeline_exp_list.txt"
awk -F'\t' '{print $1"\t"$3}' "$PARSED/COMBOS.tsv" > "$EXPLIST"
ALLMAPNAMES=$(echo "$MAPNAMES" | paste -sd,)
Rscript "$SCRIPTS/11_comparePairwiseSharedPositions_uniform.R" "$WORK/" "$ALLMAPNAMES" \
  "$(basename "$EXPLIST")" "plots/pairwise_shared_positions.pdf"

echo
echo "=== done: <mapName>.txt.annotated_tables.txt written for each construct; pairwise_shared_positions.pdf written ==="
