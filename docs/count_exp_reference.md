# Real `cut_norm.csh` / `count_exp.csh` commands (reference for stage 01 `count` mode)

`scripts/01_run_cutadapt.sh count <name>` is a parameterised template of the real per-run
`.csh` job scripts below — not a new invention. This file collects the verbatim source
commands used to build and verify it, so the correspondence stays checkable.

## cutadapt call (Norm/Exp libraries)

`IntronTRIP_run2/cut_norm.csh`:

```csh
cat names.norm.txt | while read p; do echo $p; cutadapt -j 10 -n 1 --cut 4 -g AGAAGTTCAACTCGGTACCAGAGCCACCNNNNNNNNNNNNNNNNNNGGATC --rename='{id} umi={cut_prefix} adapter={adapter_name} barcode={match_sequence}' -o ./01_cutadapt/${p}.R1.fastq.gz -p ./01_cutadapt/${p}.R2.fastq.gz ${p}_R1_001.fastq.gz ${p}_R2_001.fastq.gz --discard-untrimmed > ${p}_cutadapt.log;done
```

Note `-g` here is **unanchored** (no leading `^`). Stage 01's `count` mode had briefly
regressed to an anchored `-g "^..."`, which is stricter than this real command (requires
the primer at read position 0 exactly) — fixed to match this verbatim.

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
