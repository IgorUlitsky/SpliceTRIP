#!/bin/bash
# =============================================================================
# TRIP unified pipeline — orchestrator
# =============================================================================
# Reads a config text file (see config/example_config.txt) and runs stages 01-07
# (cutadapt -> bowtie2 -> analyzeTRIPMapping -> filterMappings[+TTAA] ->
# filterContaminants -> get_genuine_barcodes -> normalizeTRIP[+oneMismatch]) in order,
# ending at combined_counts.txt.
#
# Usage
#   bash run_pipeline.sh <config.txt> [target] [threads]
#
#   target (default: combined_counts) — how far to run:
#     cutadapt          stage 01 only
#     mapping           01-02 (Map-mode samples only)
#     genuine_barcodes  01-06
#     combined_counts   01-07 (everything)
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
  cutadapt|mapping|genuine_barcodes|combined_counts) ;;
  *) echo "unknown target '$TARGET' (expected: cutadapt | mapping | genuine_barcodes | combined_counts)" >&2; exit 1;;
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

for section in GLOBAL SAMPLES MAPS COMBOS; do
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
: "${CONT_FACTOR:=30}"
: "${R1_PRIMER:?R1_PRIMER missing from ### GLOBAL}"
: "${R1_SUFFIX:?R1_SUFFIX missing from ### GLOBAL}"
: "${R2_PRIMER:?R2_PRIMER missing from ### GLOBAL}"
TTAA_BED=${TTAA_BED:-}

echo "--- GLOBAL settings ---"
echo "BARCODE_LEN=$BARCODE_LEN  CONT_FACTOR=$CONT_FACTOR  TTAA_BED=${TTAA_BED:-<off>}"
echo "BOWTIE2_INDEX=$BOWTIE2_INDEX"
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
# Stage 01 (+ barcode/count extraction) and stage 02 (Map-mode only)
# =============================================================================
echo "=== Stage 01: cutadapt ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mode group; do
  [ -z "$name" ] && continue
  echo "-- $name (mode=$mode) --"
  stageFastq "$name" "$pathPrefix" "$r1suf" "$r2suf"
  bash "$SCRIPTS/01_run_cutadapt.sh" "$mode" "$name" "$THREADS"
done < "$PARSED/SAMPLES.tsv"
echo

if [ "$TARGET" = "cutadapt" ]; then echo "target=cutadapt reached, stopping."; exit 0; fi

echo "=== Stage 02: bowtie2 (map-mode samples only) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mode group; do
  [ -z "$name" ] && continue
  [ "$mode" = "map" ] || continue
  echo "-- $name --"
  bash "$SCRIPTS/02_run_bowtie2.sh" "$name" "$THREADS" "$BOWTIE2_INDEX"
done < "$PARSED/SAMPLES.tsv"
echo

if [ "$TARGET" = "mapping" ]; then echo "target=mapping reached, stopping."; exit 0; fi

# =============================================================================
# Stage 03/04: per map-mode sample
# =============================================================================
echo "=== Stage 03: analyzeTRIPMapping ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mode group; do
  [ -z "$name" ] && continue
  [ "$mode" = "map" ] || continue
  echo "-- $name --"
  Rscript "$SCRIPTS/03_analyzeTRIPMapping.R" "$WORK/" "$WORK/02_mapping/" "$name" "$BARCODE_LEN"
done < "$PARSED/SAMPLES.tsv"
echo

echo "=== Stage 04: filterMappings (+ TTAA if configured) ==="
while IFS=$'\t' read -r name pathPrefix r1suf r2suf mode group; do
  [ -z "$name" ] && continue
  [ "$mode" = "map" ] || continue
  echo "-- $name --"
  if [ -n "$TTAA_BED" ]; then
    Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$WORK/" "$name" ttaaFilter "$TTAA_BED"
  else
    Rscript "$SCRIPTS/04_filterMappings_TTAA.R" "$WORK/" "$name"
  fi
done < "$PARSED/SAMPLES.tsv"
echo

# =============================================================================
# Stage 05: filterContaminants -- one run PER contamination batch (7th SAMPLES column).
# Index-hopping contamination is a sequencing-run/flowcell artefact: samples that were
# never multiplexed together can't contaminate each other, so real historical fileLists
# (e.g. run3.txt/run5.txt, IntronTRIP_run2.txt/_run3.txt/_run4.txt) are scoped one per
# physical run, not one global fileList across every map sample in a combined dataset.
# Samples with no batch value fall into a single "default" batch (old/simple configs
# without a 7th column keep working as one global run, same as before this change).
# =============================================================================
echo "=== Stage 05: filterContaminants ==="
BATCHES=$(awk -F'\t' '$5=="map"{print ($7==""?"default":$7)}' "$PARSED/SAMPLES.tsv" | sort -u)
if [ -z "$BATCHES" ]; then
  echo "no map-mode samples, skipping"
else
  while read -r batch; do
    [ -z "$batch" ] && continue
    FILELIST="$WORK/.pipeline_fileList.${batch}.txt"
    awk -F'\t' -v b="$batch" '$5=="map" && ($7==""?"default":$7)==b {print $1"\t"$6}' "$PARSED/SAMPLES.tsv" > "$FILELIST"
    echo "-- batch: $batch ($(wc -l < "$FILELIST") sample(s)) --"
    Rscript "$SCRIPTS/05_filterContaminants_uniform.R" "$WORK/" "$(basename "$FILELIST")" "$CONT_FACTOR"
  done <<< "$BATCHES"
fi
echo

# =============================================================================
# Stage 06: get_genuine_barcodes -- one run per ### MAPS row (one-to-many manifest)
# =============================================================================
echo "=== Stage 06: get_genuine_barcodes ==="
while IFS=$'\t' read -r mapName sampleList; do
  [ -z "$mapName" ] && continue
  manifest="$WORK/${mapName}.txt"
  : > "$manifest"
  IFS=',' read -ra samples <<< "$sampleList"
  for s in "${samples[@]}"; do
    echo "${s}.retained_barcodes.bed" >> "$manifest"
  done
  echo "-- $mapName (${#samples[@]} sample(s)) --"
  Rscript "$SCRIPTS/06_get_genuine_barcodes.R" "$WORK/" "$(basename "$manifest")" "$BARCODE_LEN"
done < "$PARSED/MAPS.tsv"
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
echo "=== done: combined_counts.txt written for each ### COMBOS row ==="
