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
`<name>.barcodes_poses.txt` (position + strand + `barcode:count,...` list +
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
`<name>.filtered_barcodes_poses.txt` (11 columns — stage 03's 9 plus `MaxCount`/`TopBarcode`)
and `<name>.filtered_barcodes.bed`.

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
count is within `contFactor` (default 10 via the unified pipeline's `CONT_FACTOR`, i.e.
≥1/10th) of that cross-group max, or if the
key was never observed outside its own group at all (kept, since there's nothing to compare
against). In the unified pipeline, `run_pipeline.sh` builds one fileList per
`contamination_batch` value (see "Production configs" below) rather than one global list,
since contamination can only happen between samples multiplexed on the same physical
sequencing run. Output: `<name>.retained_barcodes_counts.txt` and
`<name>.retained_barcodes.bed` (the two are filtered independently, since a bed row's `V4`
field can pack several barcode variants at one position — see the script's own comments for
a real row-misalignment bug this fixes over the simpler upstream version).

### [06] get_genuine_barcodes — `scripts/06_get_genuine_barcodes.R <path> <mapName>.txt [barcodeLen] [dpniiBedPath]`

Combines every replicate Map-mode sample belonging to one construct (every `### FASTQ_MAPS`
row sharing a `mapName`, e.g. `5ss_Map_rep1`+`5ss_Map_rep2` -> mapName `5ss`) into one
genuine-barcode map for that construct. Reads each sample's `.retained_barcodes.bed`,
re-aggregates by `(Chr, Start, End, Strand)` across ALL of them, and within each position collapses near-duplicate
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

**DpnII distance.** For every genuine (singleton) position, also computes the distance to
the nearest DpnII site *downstream* of it (downstream relative to the integration's own
strand: `+` looks for the next site at a higher coordinate, `-` at a lower one), against
`DPNII_BED`/`dpniiBedPath` (default `/home/projects/ulitsky/igoru/data/human/hg38_dpnii.bed`
if not given — a Hi-C-style DpnII digest-fragment map, not a direct site list; site
positions are recovered as each fragment's start, excluding each chromosome's first
fragment, which starts at 0 as a chromosome-boundary artefact, not a real cut site). Written
to `<mapName>.txt.retained_barcodes_genuine.dnpII_dists.txt` (`chr, start, strand, barcode,
dpnII_dist`, with a header row — unlike this stage's other outputs, which are headerless
BED). Verified on a real 5ss construct: median downstream distance 179 bp, consistent with
the expected ~256 bp average genomic DpnII spacing (median of an exponential distance-to-
nearest-point-in-one-direction ≈ 0.69 × mean spacing ≈ 177 bp).

Compare against a random-genomic background with `scripts/plot_dpnII_distances.R
<path> <distFiles> <ttaaBedPath> <dpniiBedPath> <outPdf> [nSample=10000]` — samples
`nSample` random sites from the real hg38 TTAA reference (each given a random `+`/`-`
strand, since piggyBac has no inherent insertion-orientation preference), computes the same
downstream-DpnII-distance metric for each, and plots the two distributions (real vs random)
overlaid, at both linear (capped at the pooled 99th percentile) and log10 scale. `distFiles`
is a comma-separated list of `*.dnpII_dists.txt` filenames to pool into one histogram. On
the same real 5ss construct, the random-TTAA background's median was 308 bp vs the real
barcodes' 179 bp — real integration positions sit measurably closer to DpnII sites than
random genomic sites do, worth investigating further (e.g. a library-prep or mapping bias
toward reads near cut sites) rather than assuming it's expected.

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

### [08] make_tracks — `scripts/08_make_tracks_genuine_barcodes.sh <name_list.txt> <mane_bed> <mane_introns_bed> <gencode_noncoding_bed> <repeats_bed>`

Reference annotation, independent of stage 07 (only needs stage 06's genuine-barcode
maps). One `bedtools intersect` per (construct × reference track × strand orientation),
for every construct listed in `<name_list.txt>` (one `<mapName>.txt` per line — the same
format as the real production `names.complete_maps.txt`, and what the orchestrator builds
automatically from `### FASTQ_MAPS`'s distinct `mapName`s). Produces, per construct:
`<name>.MANE_all.txt`, `.MANE_exonic.txt`, `.MANE_all_antisense.txt`,
`.MANE_exonic_antisense.txt` (all against `mane_bed`, sense/antisense via `-s`/`-S`, exonic
via `-split`), `.MANE_introns_sense.txt`/`.MANE_introns_antisense.txt` (against
`mane_introns_bed`), `.gencode_noncoding_all.txt`, and `.repeats.txt` — all read by stage 10.

This is igoru's `make_tracks_genuine_barcodes.sh` verbatim, with its 4 previously-hardcoded
annotation bed paths turned into required arguments (config-driven `MANE_BED`,
`MANE_INTRONS_BED`, `GENCODE_NONCODING_BED`, `REPEATS_BED`) and `cat file | while read`
replaced with `while read < file` (same behaviour, one fewer process). Verified
byte-for-byte identical output against the original script on a real construct (5ss, 3292
genuine barcodes).

### [09] gc_content — `scripts/09_gc_content_single_genuine_barcodes.sh <name> <chrom_sizes> <genome_fasta>`

Run once **per construct** (unlike 08/10, which take the whole name-list at once) —
matches the original per-construct `bsub` usage. Slops each genuine barcode's position
±500bp (`bedtools slop`), extracts the flanking sequence (`bedtools getfasta`), computes
%GC per barcode, and writes `<name>_with_gc.bed` (stage 06's bed + a GC% column), which
stage 10 reads for `GC`/`GC_category` (Low/Average/High, ±1 SD from the construct's mean).

This is igoru's `gc_content_single_genuine_barcodes.sh`, with its 2 previously-hardcoded
reference paths turned into required arguments (config-driven `CHROM_SIZES`,
`GENOME_FASTA`) and the inline `ml BEDTools` dropped (this repo's convention: modules are
loaded once up front by the caller — see "Running" below — not re-loaded per stage script).
Verified byte-for-byte identical output against the original script.

### [10] table_annotation — `scripts/10_table_annotation_genuine_barcodes.R <path> <name_list.txt> [mane_bed] [mane_introns_bed] [gencode_noncoding_bed] [trans2gene] [mane_expression] [chromhmm_bed] [spin_bed]`

The per-construct annotation table stage 10-and-later analysis scripts (e.g.
`annotateTRIP_positions_genuine_barcodes.R`, `compareReplicatesAndMultiIntegrations_uniform.R`)
read as `<mapName>.txt.annotated_tables.txt`. For every construct in `<name_list.txt>`,
joins stage 06's genuine-barcode bed with stage 08's MANE/intron/gencode/repeats tracks and
stage 09's GC content, then adds: host-gene identity/expression (sense + antisense, exonic
+ intronic, normalized to K562 RNA-seq via `trans2gene`/`mane_expression`), distance-to-
TSS/TES/intron-boundary and relative position within the host feature, and two independent
chromatin-state annotations per position — ENCODE ChromHMM (`chromhmm_bed`) and SPIN
nuclear-compartment state (`spin_bed`).

This is igoru's `table_annotation_genuine_barcodes.R`, with its 7 previously-hardcoded
annotation paths turned into optional CLI args (each keeps its original hardcoded path as
a default when omitted, same length-guarded pattern as stage 06's `barcodeLen`/
`dpniiBedPath`) so the config drives them instead — behaviour is otherwise unchanged
(including a commented-out, never-active linear-regression residuals block from the
original, removed here as dead code rather than ported). Verified byte-for-byte identical
output against the original script on a real construct.

## Running

```bash
ml R/4.3.1-foss-2022b cutadapt/4.2-GCCcore-11.3.0 Bowtie2/2.5.1-GCC-11.3.0 \
   SAMtools/1.16.1-GCC-11.3.0 BEDTools/2.30.0-GCC-11.2.0

mkdir my_run && cd my_run
bash ../run_pipeline.sh ../config/example_config.txt combined_counts
```

`target` (2nd arg) controls how far to run: `cutadapt` | `mapping` | `genuine_barcodes` |
`combined_counts` (default — stages 01-07, the pipeline's main deliverable) | `annotated`
(stages 01-10, adds per-construct `annotated_tables.txt`). Everything is written into the
current directory. This runs stages sequentially/locally; stages 01-02 are heavy on real
libraries (tens of minutes, several GB) and should be wrapped in `bsub` for production use.

See `config/example_config.txt` for the config format:

- `### GLOBAL` — `key<TAB>value` pairs. Core: `BARCODE_LEN`, `BOWTIE2_INDEX`, `CONT_FACTOR`,
  `TTAA_BED`, `R1_PRIMER`, `R1_SUFFIX`, `R2_PRIMER`. Stages 08-10 (`annotated` target only,
  all shared cluster reference resources — safe to leave unset if you only need
  `combined_counts`): `MANE_BED`, `MANE_INTRONS_BED`, `GENCODE_NONCODING_BED`,
  `REPEATS_BED`, `CHROM_SIZES`, `GENOME_FASTA`, `TRANS2GENE`, `MANE_EXPRESSION`,
  `CHROMHMM_BED`, `SPIN_BED`.
- `### FASTQ_MAPS` — one row per Map-mode (DNA) fastq sample: `name, path_prefix, r1_suffix,
  r2_suffix, mapName, contamination_batch(optional)`. `mapName` is both which construct's genuine
  map this sample feeds into (stage 06) and its stage 05 contamination group — the two are the
  same concept here.
- `### FASTQ_NORM` — one row per Norm-mode fastq sample: `name, path_prefix, r1_suffix, r2_suffix`.
- `### FASTQ_SAMPLES` — one row per Exp-mode fastq sample, same shape as `### FASTQ_NORM`; kept as
  a separate section purely so it's clear at a glance which samples are Norm vs Exp.
- `### COMBOS` — `mapName, normSample, expSample, effLength, reverseComplement,
  oneMismatch(optional)`, referencing names from the three sections above.

Verified end-to-end on real (subsampled) data from stages 01 through 07.

## Production configs

- `config/single_combined_config.txt` — 4 constructs (5ss, CPA, 5ss_noP, CPA_noP): 12 Map + 8 Norm
  + 48 Exp samples, 48 combos. Derived verbatim from `single_combined/{5ss,CPA,5ss_noP,CPA_noP}.txt`
  + `combos.txt`.
- `config/introntrip_combined_config.txt` — 8 constructs (A1, B1, C1, D1, E1, H1, B3-6, B3-9):
  32 Map + 16 Norm + 48 Exp samples, 48 combos. Derived verbatim from
  `IntronTRIP_combined/{A1,...,B3-9}.txt` + `combos.txt`.

Both configs use the optional 6th `### FASTQ_MAPS` column, `contamination_batch`, so stage 05 runs
filterContaminants once per real sequencing run (matching the historical `run3.txt`/`run5.txt` /
`IntronTRIP_run{2,3,4}.txt` fileLists) rather than pooling every map sample across runs into one
comparison — index-hopping contamination can only happen between samples multiplexed together on
the same flowcell.

See `docs/count_exp_reference.md` for the real historical `cut_norm.csh`/`count_exp.csh` commands
that stage 01's `count` mode is a parameterised template of.

Both configs also set `DPNII_BED` (see stage 06 above) to the real hg38 DpnII digest map.

## Resolved from the two repos' open discrepancies

1. **Stage 03 `barcodes_poses.txt` column count.** ~~Kept the alignment-`Length` column (10 cols
   from stage 03, 12 after stage 04 adds `MaxCount`/`TopBarcode`)~~ — **also retracted.** Length
   was originally kept to match `ArtyomTyukaev/trip-pipeline`'s version, but it isn't just an extra
   informational column: it was also added to stage 03's position-clustering `aggregate()` grouping
   (`by = chrom, Start, strand, Length` instead of `chrom, Start, strand`), which changes what counts
   as "the same position" — two reads at the same position+strand but different alignment length
   became separate groups. igoru's own `analyzeTRIPMapping_uniform.R` never had this column, and
   position+strand (no length, no barcode) is enough to define a group. Dropped entirely: stage 03 no
   longer computes or groups by `Length`, and stage 05 is back to its original 9/11-column
   expectation (`Barcodes` col 9, `MaxCount` col 10) matching `ArtyomTyukaev/trip-pipeline`'s
   unmodified stage 05 too.
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
