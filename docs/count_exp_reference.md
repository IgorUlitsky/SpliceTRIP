# Real `cut_norm.csh` / `count_exp.csh` commands (reference for stage 01 `count` mode)

`scripts/01_run_cutadapt.sh count <name>` is a parameterised template of real production
job scripts — not a new invention. This file collects the verbatim source commands used to
build and verify it, so the correspondence stays checkable.

## cutadapt call (Norm/Exp libraries): anchored `-g`

Stage 01's `count` mode anchors `-g` (`-g "^${R1_PRIMER}${NS}${R1_SUFFIX}"`), per igoru's
explicit direction. This matches the real, actually-*submitted* per-sample bsub `.csh` job
scripts, e.g.:

`run2/5ss_Norm_rep1_S60.csh`:
```csh
bsub ... "ml cutadapt; cutadapt -j 20 -n 1 --cut 4 -g ^AGAAGTTCAACTCGGTACCAGAGCCACCNNNNNNNNNNNNNNNNNN --rename='{id} umi={cut_prefix} adapter={adapter_name} barcode={match_sequence}' -o 01_cutadapt/5ss_Norm_rep1_S60.R1.fastq.gz -p 01_cutadapt/5ss_Norm_rep1_S60.R2.fastq.gz 5ss_Norm_rep1_S60_R1_001.fastq.gz 5ss_Norm_rep1_S60_R2_001.fastq.gz --discard-untrimmed > 5ss_Norm_rep1_S60_cutadapt.log"
```

`IntronTRIP_run2/TRIP_RNA_spliced_B3-6_1-1_S53.csh` (anchored, and *with* the `GGATC`
suffix that stage 01 also includes):
```csh
bsub ... "ml cutadapt; cutadapt -j 20 -n 1 --cut 4 -g ^AGAAGTTCAACTCGGTACCAGAGCCACCNNNNNNNNNNNNNNNNNNGGATC -G  ^CGGGTTACAAAGGACTCAGGATCCACGAACTCAGCAGTGATGATTATTCTACAGAAGAAGA --rename='{id} umi={cut_prefix} adapter={adapter_name} barcode={match_sequence}' -o TRIP_RNA_spliced_B3-6_1-1_S53.R1.fastq.gz -p TRIP_RNA_spliced_B3-6_1-1_S53.R2.fastq.gz TRIP_RNA_spliced_B3-6_1-1_S53_R1_001.fastq.gz TRIP_RNA_spliced_B3-6_1-1_S53_R2_001.fastq.gz --discard-untrimmed > TRIP_RNA_spliced_B3-6_1-1_S53_cutadapt.log"
```

**Real historical usage is inconsistent across runs** — the generic `cut_norm.csh`/
`count_exp.csh` batch-loop templates in some run directories (e.g. `IntronTRIP_run2/
cut_norm.csh`, below) use an *unanchored* `-g`, and a few per-sample scripts
(`IntronTRIP_run4/Norm_B1_1_1_S1.csh`, `SON_run2/SON_TRIP_A1_1_1_DMSO_Spliced_Rep2_S24.csh`)
are also unanchored. Anchored was chosen as the convention to follow going forward.

`IntronTRIP_run2/cut_norm.csh` (unanchored batch-loop template, shown for contrast — not
what stage 01 follows):

```csh
cat names.norm.txt | while read p; do echo $p; cutadapt -j 10 -n 1 --cut 4 -g AGAAGTTCAACTCGGTACCAGAGCCACCNNNNNNNNNNNNNNNNNNGGATC --rename='{id} umi={cut_prefix} adapter={adapter_name} barcode={match_sequence}' -o ./01_cutadapt/${p}.R1.fastq.gz -p ./01_cutadapt/${p}.R2.fastq.gz ${p}_R1_001.fastq.gz ${p}_R2_001.fastq.gz --discard-untrimmed > ${p}_cutadapt.log;done
```

## Barcode counting (the "`count_exp.csh`" step)

Same one-liner, reused across every run directory checked, differing only in the
`substr($4, OFFSET, 18)` offset (`OFFSET` = `len(R1_PRIMER) + 1`):

`run2/count_exp.csh` (28-nt primer, offset 29 — also processes `names.norm.txt`):
```csh
cat names.norm.txt names.exp.txt | while read p; do echo $p; zcat ./01_cutadapt/${p}.R1.fastq.gz | sed 's/umi=//' | sed 's/barcode=//' | awk '{if (NR % 4 == 1) {print substr($4,29,18)}}' | sort | uniq -c | sort -k1,1rn | awk '{print $2 "\t" $1}' | gzip > ${p}.counts.txt.gz;zcat ${p}.counts.txt.gz | wc -l; done
```

`run4/count_exp.csh` (28-nt primer, offset 29, `names.norm.txt` only):
```csh
cat names.norm.txt | while read p; do echo $p; zcat ./01_cutadapt/${p}.R1.fastq.gz | sed 's/umi=//' | sed 's/barcode=//' | awk '{if (NR % 4 == 1) {print substr($4,29,18)}}' | sort | uniq -c | sort -k1,1rn | awk '{print $2 "\t" $1}' | gzip > ${p}.counts.txt.gz;zcat ${p}.counts.txt.gz | wc -l; done
```

`run7/count_exp.csh` (**59**-nt offset — run7's RNA-readout libraries use the longer 59-nt
`R1_PRIMER`, confirming stage 01's "28 -> 29, 58 -> 59" offset comment):
```csh
cat names.exp.txt | while read p; do echo $p; zcat ./01_cutadapt/${p}.R1.fastq.gz | sed 's/umi=//' | sed 's/barcode=//' | awk '{if (NR % 4 == 1) {print substr($4,59,18)}}' | sort | uniq -c | sort -k1,1rn | awk '{print $2 "\t" $1}' | gzip > ${p}.counts.txt.gz;zcat ${p}.counts.txt.gz | wc -l; done
```

`IntronTRIP_run2/count_exp.csh` and `SON_run2/count_exp.csh`: identical to `run2`'s /
`run4`'s form (offset 29).

Stage 01's `count` mode combines the cutadapt call and this counting one-liner into a
single `bash 01_run_cutadapt.sh count <name>` invocation, with `R1_PRIMER`/`OFFSET`
parameterised instead of hard-coded per run.

## Not adopted: `IntronTRIP_run2/count_norm.csh`

This one script is a one-off variant, not the general convention — it has a commented-out
UMI-extraction line and an *active* line that assumes a pre-existing 3-column
`barcodes.txt.gz` and does a `cut -f2-3 | sort | uniq` UMI-dedup step before counting.
Every other `count_exp.csh` across `run2`/`run4`/`run7`/`SON_run2`/`IntronTRIP_run2` uses
the simple read-count form above (no UMI dedup for Norm/Exp libraries), which is what
stage 01 implements.
