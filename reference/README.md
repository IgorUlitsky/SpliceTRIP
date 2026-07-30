# Bundled reference tracks

Genome annotation tracks used by the `annotated` target (stages 08-10) that are small
enough to bundle directly in the repo, for anyone running this pipeline without access to
the Weizmann cluster's shared storage. Point the matching `### GLOBAL` config key at the
path here instead of the cluster default to use them.

| File | Config key | What it is | Source |
|---|---|---|---|
| `MANE.GRCh38.v1.0.refseq_genomic.bed` | `MANE_BED` | MANE (Matched Annotation from NCBI and EMBL-EBI) v1.0 RefSeq gene models, hg38, BED format | [NCBI MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/), release v1.0 |
| `K562_ChromHMM_Kellis_ENCFF319VXX.bed.gz` | `CHROMHMM_BED` | ChromHMM chromatin-state segmentation for K562 (Kellis lab), hg38, gzipped BED | [ENCODE](https://www.encodeproject.org/), accession `ENCFF319VXX` |
| `K562_SPIN_25kb_hg38.bed` | `SPIN_BED` | SPIN (Spatial Position Inference of the Nuclear genome) nuclear-compartment calls for K562, 25kb bins, hg38 | Spin-state track, see [Wang et al. 2021](https://doi.org/10.1101/gr.266239.120) |

The remaining `annotated`-target reference keys (`GENOME_FASTA`, `CHROM_SIZES`,
`MANE_INTRONS_BED`, `GENCODE_NONCODING_BED`, `REPEATS_BED`, `TRANS2GENE`,
`MANE_EXPRESSION`) are either too large (the genome FASTA) or too environment-specific to
bundle here — point them at your own copies of the equivalent resources.
