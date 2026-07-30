# Sample `expListFile`s for stage 11 (comparePairwiseSharedPositions)

Real, verbatim production `expListFile`s (as consumed by
`scripts/11_comparePairwiseSharedPositions_uniform.R`'s 3rd argument), included here as
concrete examples of the format and of real construct/expression-library groupings.

Format: `<mapName><TAB><expName>`, one row per (construct, expression library) pair.
`<expName>` must have a matching `<expName>.combined_counts.txt` next to the `<mapName>.txt.
annotated_tables.txt` files (i.e. produced by stage 07, or by an equivalent pipeline that
writes the same columns — see below).

- **`single_combined.exp_list.txt`** — 4 constructs (`5ss`, `CPA`, `5ss_noP`, `CPA_noP`), 2
  expression replicates each. `exp.log10` here is the standard stage 07
  `log10(1+exp)` — these `.combined_counts.txt` come straight from this pipeline's own
  stage 07.

- **`IntronTRIP_combined.spliced_unspliced.exp_list.txt`** and
  **`SON_combined.spliced_unspliced.exp_list.txt`** — every `<expName>` ends in
  `.spliced_unspliced`. These `.combined_counts.txt` do **not** come from stage 07 — they're
  produced by a separate splicing-ratio pipeline where `exp.log10` is a bounded [0,1]
  spliced/(spliced+unspliced) ratio, not `log10(1+exp)`. Stage 11 detects this
  automatically (it only computes `log10(1+exp)` when the file doesn't already have its own
  `exp.log10` column) and uses the ratio as-is — this is exactly the scenario that fix was
  written for (see stage 11's script comments, and the git log for the bug this caught:
  values >1 in a column that should never exceed 1).

## Example invocation

```bash
Rscript scripts/11_comparePairwiseSharedPositions_uniform.R \
  /path/to/IntronTRIP_combined/ \
  A1,B1,C1,D1,E1,B3-6,B3-9,H1 \
  examples/exp_lists/IntronTRIP_combined.spliced_unspliced.exp_list.txt \
  plots/comparePairwiseSharedPositions_IntronTRIP_splicing.pdf \
  splicing
```

(`examples/exp_lists/...` only resolves correctly if you copy/symlink the file into
`/path/to/IntronTRIP_combined/` first — stage 11 reads `<path><expListFile>`, same
convention as every other stage.)
