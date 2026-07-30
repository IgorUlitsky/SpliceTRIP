#!/bin/bash
# =============================================================================
# Stage 17 — metagene bigWig overlap prep (histone marks + speckle/TSA-seq tracks)
# =============================================================================
# Sequential port of the real meta_histone_all.csh + meta_histone_summary.sh +
# meta_speckle_all.csh + meta_speckle_summary.sh (single_combined/*.csh/*.sh) -- same
# commands, just run directly one after another instead of each bigWigAverageOverBed call
# being its own async bsub job with LSF job-ID dependency tracking (this pipeline's own
# convention: it doesn't submit jobs itself, see run_pipeline.sh's own header comment).
# For production-scale runs (many histone marks/speckle tracks x many constructs), wrap
# this whole stage in bsub yourself rather than running it on a login node.
#
# Usage
#   bash 17_meta_bigwig_overlap.sh <mapTypesFile> <histoneBwDir> <speckleDir> <chromSizes> <chainHg38ToHg19>
#
# Requires bedtools + kentUtils modules loaded (bigWigAverageOverBed, liftOver).
#
# Inputs   <mapTypesFile>              one bare mapName per line (matching ### FASTQ_MAPS
#                                      mapNames -- run_pipeline.sh generates this)
#          <mapName>.txt.retained_barcodes_genuine.bed         (stage 06)
#          histone_marks.txt           one histone/factor track name per line -- a real,
#                                      curated per-dataset file (see e.g.
#                                      single_combined/histone_marks.txt), not generated
#                                      here; must already exist in the working directory
#          speckle_names.txt           same, for speckle/TSA-seq tracks
#          <histoneBwDir>/<mark>.K562.bw       ENCODE K562 histone/factor bigWigs (hg38)
#          <speckleDir>/<track>.bw             speckle/TSA-seq bigWigs (hg19)
#
# Outputs  <mapName>.plus_50_each_side.bed          bedtools slop -b 50 window per barcode
#          hg19.<mapName>.plus_50_each_side.bed     liftOver'd version, for the hg19 speckle tracks
#          meta/<mapName>.<mark>.tab.gz             raw bigWigAverageOverBed output, one per (construct, track)
#          meta/<mapName>.combined_histone.txt      matrix: one row per barcode, one column per
#                                                    histone_marks.txt entry (alphabetical), value = mean0
#          meta/<mapName>.combined_speckle.txt      same, for speckle_names.txt
#          (what scripts/18_plotMetageneCorrelations_uniform.R reads)
# =============================================================================
set -euo pipefail

MAPTYPES_FILE=${1:?"mapTypesFile required (one mapName per line)"}
HISTONE_BW_DIR=${2:-/home/projects/ulitsky/shared/data/human/ENCODE/K562}
SPECKLE_DIR=${3:-/home/projects/ulitsky/shared/data/human/speckles}
CHROM_SIZES=${4:-/home/projects/ulitsky/shared/ucsctracks/hg38/hg38.chrom.sizes}
CHAIN=${5:-/home/projects/ulitsky/igoru/data/combinedTables/chains/hg38ToHg19.over.chain}

mkdir -p meta

echo "-- bedtools slop -b 50 per construct --"
while read -r p; do
  [ -z "$p" ] && continue
  bedtools slop -b 50 -i "${p}.txt.retained_barcodes_genuine.bed" -g "$CHROM_SIZES" > "${p}.plus_50_each_side.bed"
done < "$MAPTYPES_FILE"

echo "-- histone/factor tracks: bigWigAverageOverBed (hg38) --"
while read -r p; do
  [ -z "$p" ] && continue
  while read -r q; do
    [ -z "$q" ] && continue
    bigWigAverageOverBed "${HISTONE_BW_DIR}/${q}.K562.bw" "${p}.plus_50_each_side.bed" "meta/${p}.${q}.tab"
    gzip -f "meta/${p}.${q}.tab"
  done < histone_marks.txt
done < "$MAPTYPES_FILE"

echo "-- liftOver to hg19 (speckle/TSA-seq tracks are hg19) --"
while read -r p; do
  [ -z "$p" ] && continue
  liftOver "${p}.plus_50_each_side.bed" "$CHAIN" "hg19.${p}.plus_50_each_side.bed" unmapped.txt
done < "$MAPTYPES_FILE"

echo "-- speckle/TSA-seq tracks: bigWigAverageOverBed (hg19) --"
while read -r p; do
  [ -z "$p" ] && continue
  while read -r q; do
    [ -z "$q" ] && continue
    bigWigAverageOverBed "${SPECKLE_DIR}/${q}.bw" "hg19.${p}.plus_50_each_side.bed" "meta/${p}.${q}.tab"
    gzip -f "meta/${p}.${q}.tab"
  done < speckle_names.txt
done < "$MAPTYPES_FILE"

# Column 5 of bigWigAverageOverBed's output is mean0 (mean over the region, treating
# uncovered bases as 0) -- confirmed against real production meta/*.tab.gz files.
# Globbing meta/$p.*.column_*.txt for paste (not iterating histone_marks.txt/
# speckle_names.txt again) so the column order is alphabetical by track name, matching
# the R script's own sort(histone_marks.txt)/sort(speckle_names.txt).
echo "-- combining into meta/<mapName>.combined_{histone,speckle}.txt --"
while read -r p; do
  [ -z "$p" ] && continue
  while read -r q; do zcat "meta/${p}.${q}.tab.gz" | cut -f5 > "meta/${p}.${q}.column_histone.txt"; done < histone_marks.txt
  paste meta/${p}.*.column_histone.txt > "meta/${p}.combined_histone.txt"
  while read -r q; do zcat "meta/${p}.${q}.tab.gz" | cut -f5 > "meta/${p}.${q}.column_speckle.txt"; done < speckle_names.txt
  paste meta/${p}.*.column_speckle.txt > "meta/${p}.combined_speckle.txt"
done < "$MAPTYPES_FILE"

echo "stage 17 done"
