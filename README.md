# trip-unified-pipeline

Config-driven TRIP primary-processing pipeline: raw paired fastq -> `combined_counts.txt`.

Built by merging igoru's own scripts with the partially-revised versions from
[ArtyomTyukaev/trip-pipeline](https://github.com/ArtyomTyukaev/trip-pipeline), resolving the two
open discrepancies between them (see below), and adding a single config file + orchestrator on top
so a run is driven by one text file instead of hand-written per-sample shell/bsub commands.

## Pipeline

```
raw fastq (R1/R2)
 -> [01] cutadapt (two-pass: barcode/UMI into read header, then read-through trim)
 -> [02] bowtie2 mapping -> second-read BED12                    (map-mode samples only)
 -> [03] analyzeTRIPMapping   -> <name>.barcodes_poses.txt / .barcodes.bed
 -> [04] filterMappings (+ optional TTAA filter) -> <name>.filtered_barcodes_poses.txt / .bed
 -> [05] filterContaminants   -> <name>.retained_barcodes_counts.txt / .bed
 -> [06] get_genuine_barcodes -> <mapName>.txt.retained_barcodes_genuine.bed   (one-to-many: combines replicate samples)
 -> [07] normalizeTRIP (+ optional 1-nt mismatch tolerance) -> <expName>.combined_counts.txt  ★
```

## Running

```bash
ml R/4.3.1-foss-2022b cutadapt/4.2-GCCcore-11.3.0 Bowtie2/2.5.1-GCC-11.3.0 \
   SAMtools/1.16.1-GCC-11.3.0 BEDTools/2.30.0-GCC-11.2.0

mkdir my_run && cd my_run
bash ../run_pipeline.sh ../config/example_config.txt combined_counts
```

`target` (2nd arg) controls how far to run: `cutadapt` | `mapping` | `genuine_barcodes` |
`combined_counts` (default, runs everything). Everything is written into the current directory.
This runs stages sequentially/locally; stages 01-02 are heavy on real libraries (tens of minutes,
several GB) and should be wrapped in `bsub` for production use.

See `config/example_config.txt` for the config format (`### GLOBAL` / `### SAMPLES` / `### MAPS` /
`### COMBOS`). Verified end-to-end on real (subsampled) data from stages 01 through 07.

## Resolved from the two repos' open discrepancies

1. **Stage 03 `barcodes_poses.txt` column count.** Kept the alignment-`Length` column (10 cols from
   stage 03, 12 after stage 04 adds `MaxCount`/`TopBarcode`) rather than dropping it. Stage 05's
   `ncol` check and column-index references (`Barcodes` col 9->10, `MaxCount` col 10->11) were
   updated to match — a plain "just check for 12" fix alone would have silently read the wrong
   columns.
2. **`analyzeTRIPMapping` minus-strand position: "V3-1", not "V3".** `Start = V3-1` (not `V3`) for
   `-`-strand alignments — BED end is exclusive, so `V3-1` is the actual last aligned base.
   Verified empirically against the real hg38 TTAA reference on real data: with `V3-1`, stage 04's
   minus-strand TTAA offset must be `tt = Start+1` (not `tt = Start`), or the TTAA filter drops
   ~99.9% of genuine minus-strand integrations instead of the expected ~94% hit rate.

## Notable differences from `ArtyomTyukaev/trip-pipeline`

- Stage 05 here is built on igoru's current (more advanced) `filterContaminants_uniform.R`, which
  additionally fixes a real `.filtered_barcodes.bed` row-misalignment bug that Artyom's simpler copy
  doesn't have (see the script's own comments).
- Stage 05 now takes an optional 3rd CLI arg for `contFactor` so the config can drive it.
- Stage 06 still uses the older, non-vectorized `get_genuine_barcodes_v2.R` algorithm (matching
  Artyom's choice) rather than the more-refactored `get_genuine_barcodes_IntronTRIP.R` variant — a
  candidate follow-up if stage 06 becomes a bottleneck at real scale.
