# Sample `expListFile`s for comparePairwiseSharedPositions and compareReplicatesAndMultiIntegrations

Real, verbatim production `expListFile`s (as consumed by
`scripts/comparePairwiseSharedPositions_uniform.R`'s and
`scripts/compareReplicatesAndMultiIntegrations_uniform.R`'s 3rd argument — both scripts
use the identical `<mapName><TAB><expName>` format), included here as concrete examples of
the format and of real construct/expression-library groupings.

Format: `<mapName><TAB><expName>`, one row per (construct, expression library) pair.
`<expName>` must have a matching `<expName>.combined_counts.txt` next to the `<mapName>.txt.
annotated_tables.txt`/`<mapName>.txt.annotated_table_used.txt` files (i.e. produced by stage
07, or by an equivalent pipeline that writes the same columns — see below).

- **`single_combined.exp_list.txt`** — 4 constructs (`5ss`, `CPA`, `5ss_noP`, `CPA_noP`), 2
  expression replicates each. `exp.log10` here is the standard stage 07
  `log10(1+exp)` — these `.combined_counts.txt` come straight from this pipeline's own
  stage 07.

- **`IntronTRIP_combined.spliced_unspliced.exp_list.txt`** and
  **`SON_combined.spliced_unspliced.exp_list.txt`** — every `<expName>` ends in
  `.spliced_unspliced`. These `.combined_counts.txt` do **not** come from stage 07 — they're
  produced by a separate splicing-ratio pipeline where `exp.log10` is a bounded [0,1]
  spliced/(spliced+unspliced) ratio, not `log10(1+exp)`. comparePairwiseSharedPositions_uniform.R
  detects this automatically (it only computes `log10(1+exp)` when the file doesn't already
  have its own `exp.log10` column) and uses the ratio as-is — this is exactly the scenario
  that fix was written for (see the script's own comments, and the git log for the bug this
  caught: values >1 in a column that should never exceed 1).

- **`IntronTRIP_combined.expression.exp_list.txt`** — one `Total_*` readout per replicate,
  per construct (all 8 IntronTRIP_combined constructs), i.e. the real `expression.exp_list.txt`
  used with compareReplicatesAndMultiIntegrations_uniform.R for the "how well do independent
  barcodes in the same gene/intron/position agree" QC report — deliberately just the Total
  readout, not every ### COMBOS expName for these constructs (spliced/unspliced/3intron/etc
  are a different question).

## Example invocations

```bash
Rscript scripts/comparePairwiseSharedPositions_uniform.R \
  /path/to/IntronTRIP_combined/ \
  A1,B1,C1,D1,E1,B3-6,B3-9,H1 \
  examples/exp_lists/IntronTRIP_combined.spliced_unspliced.exp_list.txt \
  plots/comparePairwiseSharedPositions_IntronTRIP_splicing.pdf \
  splicing

Rscript scripts/compareReplicatesAndMultiIntegrations_uniform.R \
  /path/to/IntronTRIP_combined/ \
  A1,B1,C1,D1,E1,B3-6,B3-9,H1 \
  examples/exp_lists/IntronTRIP_combined.expression.exp_list.txt \
  plots/replicates_integrations_IntronTRIP_expression.pdf
```

(`examples/exp_lists/...` only resolves correctly if you copy/symlink the file into
`/path/to/IntronTRIP_combined/` first — both scripts read `<path><expListFile>`, same
convention as every other stage/script.)
