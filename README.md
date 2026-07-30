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

Background: a TRIP experiment transfects a pool of plasmids, each carrying a random
16- or 18-nt barcode next to a piggyBac reporter, so barcodes land at essentially random
genomic positions (`Map` libraries). The reporter's activity at each position is then read
out by sequencing the barcode alone from RNA or normalized DNA (`Norm`/`Exp` libraries).
Stages 01-06 are all about pinning down, per barcode, *where* it landed and *whether that
call can be trusted*; stage 07 is where a `Norm`/`Exp` pair actually gets turned into a
per-position expression value.

Barcode length is a real per-construct choice, not a fixed constant: `BARCODE_LEN` in
`### GLOBAL` (default 18) flows through to stage 01 (`--cut`/offset math), stage 03 (barcode
length filter), and stage 06 (barcode length filter) — set it to 16 for constructs that use
a 16-nt barcode. `effLength` in `### COMBOS` (stage 07) is a separate, per-combo setting for
how many of those bases the Norm/Exp readout actually uses, and can differ from
`BARCODE_LEN` itself.

### [01] cutadapt — `scripts/01_run_cutadapt.sh {map|count} <name>`

Extracts the UMI + barcode from R1's header into the read name (via cutadapt's `--rename`)
and trims the read-through/adapter sequence, so every later stage works off a clean,
position-tagged FASTQ instead of raw sequence.

- **`map` mode** (Map libraries — DNA, paired-end, R2 is genomic): two cutadapt passes.
  Pass 1 finds the R1 primer + `BARCODE_LEN` N's (barcode) + junction sequence with
  `--cut 4` first (strips the 4-nt UMI) and anchors R2 to the ITR primer that ends in the
  piggyBac `TTAA` duplication signature. Pass 2 trims any read-through into the construct
  from both mates. Output: `01_cutadapt/<name>.R{1,2}.fastq.gz` (R1 header carries
  `umi=`/`barcode=`) plus `<name>.barcodes.txt.gz` (read_name, UMI, barcode — one
  exactly-`BARCODE_LEN`-nt barcode per read).
- **`count` mode** (Norm/Exp libraries — barcode readout only, R2 unused downstream): a
  single **anchored** (`-g "^..."`) cutadapt pass, then an awk one-liner tallies
  `barcode -> read_count` straight from the renamed FASTQ header, sorted by count
  descending. Output: `01_cutadapt/<name>.R{1,2}.fastq.gz` and `<name>.counts.txt.gz`.
  See `docs/count_exp_reference.md` for the real per-sample job scripts this mode is a
  parameterised template of, and why the primer match is anchored.

### [02] bowtie2 mapping — `scripts/02_run_bowtie2.sh <name>` (Map-mode only)

Aligns both mates end-to-end to hg38 (`bowtie2 -X 2000`), then keeps only read 2 (flag
`0x80`) — R1 is the constant reporter sequence whose barcode was already captured in stage
01, R2 is the genomic DNA flanking the integration site. Converts R2's alignments to one
BED12 line each, piping `samtools view -f 0x80` straight into `bedtools bamtobed -i stdin`
rather than writing an intermediate read-2-only BAM (nothing downstream reads a BAM, only
the BED). Output: `<name>.bam` (both mates, sorted+indexed) and
`02_mapping/<name>.second_read.bed.gz` (what stage 03 reads). This script is a
reconstruction — no per-sample script for this step survives in the run directories, so the
alignment was historically run ad hoc; treat the bowtie2 flags as authoritative, the
packaging as new.

### [03] analyzeTRIPMapping — `scripts/03_analyzeTRIPMapping.R <path> <path_to_mapping> <name> [barcodeLen]`

The core per-read-to-per-position collapse. Joins stage 01's barcode table to stage 02's
alignment BED by read name, drops alignments outside a sane length range (20-300 bp) and
any landing in the globin locus (a common plasmid-integration artefact region on chr11),
then groups reads by genomic position + strand and asks, per position: which barcode(s)
were actually seen there, and in what proportion? A position can have several barcode
"clusters" if some are just sequencing/PCR noise around the dominant one — the position is
kept only if noise doesn't swamp the signal (`screenFunc`: drop positions where too many
distinct barcodes each carry ≥1% of the max count). Two outputs, row-aligned to each other:
`<name>.barcodes_poses.txt` (position + strand + length + `barcode:count,...` list +
aggregate counts) and `<name>.barcodes.bed` (same positions in BED6 for genome browsers).
`barcodeLen` (4th arg, default 18 if omitted; `run_pipeline.sh` always passes `BARCODE_LEN`)
controls the exact-length filter applied to barcodes before they're kept.

Minus-strand position uses plain `Start = V3` (matching both igoru's own current
`analyzeTRIPMapping_uniform.R` and the original, unmodified `trip-pipeline` script). An
earlier version of this repo used `V3-1` instead, on the theory that BED end is exclusive
so `V3-1` is the "true last aligned base" — reverted after real-data TTAA verification
showed it was wrong; see stage 04 and "Resolved discrepancies" below.

### [04] filterMappings (+ optional TTAA filter) — `scripts/04_filterMappings_TTAA.R <path> <name> [ttaaFilter <ttaa_bed>]`

Resolves multiple nearby integration calls that are really the same physical integration
site read slightly differently. Per chromosome, positions are visited in order of
decreasing read support (`MaxCount`/`TopBarcode`, computed here from stage 03's barcode
list); a position is dropped as redundant if a stronger position is within 10 bp, has
≥4x its read count, and its dominant barcode is within edit distance 3 (a near-duplicate,
not a coincidence). Survivors are written out with the two new columns appended:
`<name>.filtered_barcodes_poses.txt` (12 columns — the alignment-`Length` column from
stage 03 plus `MaxCount`/`TopBarcode`) and `<name>.filtered_barcodes.bed`.

Optional **TTAA filter** (piggyBac integrates exclusively at `TTAA` sites, so a position
that isn't adjacent to a real genomic `TTAA` is very likely a mapping artefact): pass
`ttaaFilter <ttaa_bed_path>` (a 6-col BED like `hg38_ttaa.bed`) as extra args to require
each surviving position to sit exactly at a known TTAA (strand-aware offset: `+` strand
checks `V2-4`, `-` strand checks `V2`, matching stage 03's plain-`V3` convention). Off by
default; enabled pipeline-wide by setting `TTAA_BED` in `### GLOBAL`.

### [05] filterContaminants — `scripts/05_filterContaminants_uniform.R <path> <fileList> [contFactor]`

Removes index-hopping contamination between mapping libraries that were sequenced
together. Takes a fileList (`sample_name<TAB>group`, one row per Map-mode sample) and, per
`(Chr, Pos, Barcode)` key, compares each sample's read count against the *maximum* count
seen for that same key in any *other* group on the same list. A key survives if its own
count is within `contFactor` (default 30, i.e. ≥1/30th) of that cross-group max, or if the
key was never observed outside its own group at all (kept, since there's nothing to compare
against). In the unified pipeline, `run_pipeline.sh` builds one fileList per
`contamination_batch` value (see "Production configs" below) rather than one global list,
since contamination can only happen between samples multiplexed on the same physical
sequencing run. Output: `<name>.retained_barcodes_counts.txt` and
`<name>.retained_barcodes.bed` (the two are filtered independently, since a bed row's `V4`
field can pack several barcode variants at one position — see the script's own comments for
a real row-misalignment bug this fixes over the simpler upstream version).

### [06] get_genuine_barcodes — `scripts/06_get_genuine_barcodes.R <path> <mapName>.txt [barcodeLen]`

Combines every replicate Map-mode sample belonging to one construct (the `### MAPS` list,
e.g. `5ss_Map_rep1`+`5ss_Map_rep2` -> mapName `5ss`) into one genuine-barcode map for that
construct. Reads each sample's `.retained_barcodes.bed`, re-aggregates by `(Chr, Start,
End, Strand)` across ALL of them, and within each position collapses near-duplicate
barcodes (edit distance ≤2) into their dominant variant, summing counts. A barcode/position
pair is "genuine" if that barcode maps to exactly one position across the whole construct
(`retained_barcodes_genuine.bed`); barcodes seen at more than one position are split off
separately (`retained_barcodes_multimappers.bed`) rather than being silently kept or
dropped. `_forUCSC` variants are the same tables with scores capped at 1000 for genome
browser display. This is what stage 07 reads as its map of barcode -> genomic position.

This is igoru's `~/R/trip/IntronTrip/get_genuine_barcodes_IntronTRIP.R`, not the older
`get_genuine_barcodes_v2.R` this repo started with — same grouping/merge logic (including
its documented quirk where a barcode merged into an already-dropped one can lose its count;
kept intentionally, not changed here), but the per-position dedup is rewritten with
`tapply`/vectors instead of `aggregate()`-on-a-dataframe-per-position, which mattered at
real scale (~140k positions, profiled bottleneck). It also guards the optional 3rd
(`barcodeLen`) argument — the old version's unconditional `as.integer(myArgs[3])` silently
became `NA` when called with only 2 args, which silently zeroed out every position's output
with no error at all.

### [07] normalizeTRIP (+ optional 1-nt mismatch tolerance) — `scripts/07_normalizeTRIP_MM.R <path> <path_norm> <mapName>.txt <normSample> <expSample> <effLength> <reverseComplement> [oneMismatch]`

The actual expression readout, run once per `### COMBOS` row. Loads the construct's genuine
barcode map (stage 06), the Norm library's barcode counts (stage 01 `count` mode, floored
at `minNormReads = 2`), and the Exp library's barcode counts, all truncated/aggregated to
the first `effLength` bases (some readouts use a shorter effective barcode than the
construct's full `BARCODE_LEN`). Keeps only barcodes that (a) map to exactly one genomic position and (b) are present
in both Norm and Exp, then writes one row per barcode: raw `norm`/`exp` counts plus the
mapped `Chr`/`Pos`/`Strand`. `reverseComplement` (`TRUE`/`FALSE`, per-combo in the config)
handles Exp libraries sequenced from the opposite strand of the barcode (e.g. some RNA
readouts). Output: `<expName>.combined_counts.txt` ★ (the pipeline's end product — per-
barcode norm/exp counts + position, ready for a `exp.norm = log2((1+exp)/(1+norm))`-style
downstream analysis) and `<expName>.bed`.

Optional **1-nt mismatch tolerance** (`oneMismatch` as the 8th arg, or a `yes` in the
config's `### COMBOS` 6th column): before intersecting with the genuine map, Norm/Exp
barcodes that are exactly one substitution away from a single genuine barcode (and not
ambiguous between two) are "snapped" onto it and their counts summed in — rescues reads
with a single sequencing error that would otherwise be silently dropped for not matching
any genuine barcode exactly. Off by default; exact-match behaviour is unchanged when the
flag is omitted.

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

## Production configs

- `config/single_combined_config.txt` — 4 constructs (5ss, CPA, 5ss_noP, CPA_noP), 68 samples,
  48 combos. Derived verbatim from `single_combined/{5ss,CPA,5ss_noP,CPA_noP}.txt` + `combos.txt`.
- `config/introntrip_combined_config.txt` — 8 constructs (A1, B1, C1, D1, E1, H1, B3-6, B3-9),
  96 samples, 48 combos. Derived verbatim from `IntronTRIP_combined/{A1,...,B3-9}.txt` + `combos.txt`.

Both configs use the optional 7th `### SAMPLES` column, `contamination_batch`, so stage 05 runs
filterContaminants once per real sequencing run (matching the historical `run3.txt`/`run5.txt` /
`IntronTRIP_run{2,3,4}.txt` fileLists) rather than pooling every map sample across runs into one
comparison — index-hopping contamination can only happen between samples multiplexed together on
the same flowcell.

See `docs/count_exp_reference.md` for the real historical `cut_norm.csh`/`count_exp.csh` commands
that stage 01's `count` mode is a parameterised template of.

## Resolved from the two repos' open discrepancies

1. **Stage 03 `barcodes_poses.txt` column count.** Kept the alignment-`Length` column (10 cols from
   stage 03, 12 after stage 04 adds `MaxCount`/`TopBarcode`) rather than dropping it, matching
   `ArtyomTyukaev/trip-pipeline`'s version. Stage 05's `ncol` check and column-index references
   (`Barcodes` col 9->10, `MaxCount` col 10->11) were updated to match — a plain "just check for 12"
   fix alone would have silently read the wrong columns. Note: igoru's own current
   `analyzeTRIPMapping_uniform.R` doesn't have this column at all (9 cols, no `Length`) — this was a
   deliberate choice to follow Artyom's addition rather than igoru's leaner original, since the extra
   column is purely informational and doesn't change any filtering logic.
2. **`analyzeTRIPMapping` minus-strand position.** ~~An earlier version of this repo used
   `Start = V3-1` instead of `V3` for `-`-strand alignments~~ — **retracted.** That change (and the
   matching `tt = Start+1` in stage 04) was verified only by internal self-consistency, not against
   real production data, and turned out to be wrong: cross-checked directly against real,
   already-computed `.barcodes.bed` files (`run3/5ss_Map_rep1_S81`, `run3/CPA_Map_rep3_S80`) and the
   real hg38 TTAA reference, plain `V3`/`tt = Start` gives a ~58-70% minus-strand TTAA hit rate
   (matching the `+`-strand baseline), while `V3-1`/`tt = Start+1` collapses to ~1.6-1.9%. Both
   igoru's own `analyzeTRIPMapping_uniform.R` and the original, unmodified `trip-pipeline` script use
   plain `V3` — there was never an actual discrepancy here to resolve.

## Notable differences from `ArtyomTyukaev/trip-pipeline`

- Stage 05 here is built on igoru's current (more advanced) `filterContaminants_uniform.R`, which
  additionally fixes a real `.filtered_barcodes.bed` row-misalignment bug that Artyom's simpler copy
  doesn't have (see the script's own comments).
- Stage 05 now takes an optional 3rd CLI arg for `contFactor` so the config can drive it.
- Stage 06 uses igoru's `get_genuine_barcodes_IntronTRIP.R` (vectorized, faster at real scale,
  with a fixed argument-parsing bug) rather than Artyom's `get_genuine_barcodes_v2.R` copy — see
  its README section above.
