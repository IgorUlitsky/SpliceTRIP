# =============================================================================
# plotSplicingEfficiencies — per-construct splicing-efficiency violin plot
# =============================================================================
# Visualization script, not a pipeline stage -- run standalone against already-produced
# stage 07/10/14 output, no run_pipeline.sh wiring or numeric prefix.
#
# Usage
#   Rscript plotSplicingEfficiencies_uniform.R <path> <expListFile> <outFile.pdf>
#
# expListFile: expName<TAB>niceName<TAB>nameMap.txt -- one row per violin, in the exact
#   order they should appear on the x-axis (no automatic reordering/grouping, unlike the
#   original script's grep-based construct sort -- put the rows in your file in display
#   order). Two rows sharing the same niceName are pooled into one violin.
#     expName  -- <path><expName>.combined_counts.txt must exist and already carry its own
#                 exp.log10 (a [0,1]-bounded splicing ratio, e.g. spliced/(spliced+unspliced)
#                 -- stage 14's *.spliced_unspliced.combined_counts.txt is the usual source).
#     niceName -- display label for this violin.
#     nameMap  -- <path><nameMap>.annotated_tables.txt must exist (stage 10 output);
#                 include the trailing ".txt" (i.e. pass e.g. "A1.txt", matching every other
#                 script in this pipeline's <mapName>.txt convention), giving
#                 "<nameMap>.annotated_tables.txt" -- e.g. nameMap="A1.txt" reads
#                 "A1.txt.annotated_tables.txt". Repeated nameMaps across rows are read once
#                 and reused.
#
# Output  <outFile.pdf> -- one violin per expListFile row (or per pooled niceName), showing
#   the distribution of per-position splicing efficiency (exp.log10, restricted to [0,1]).
#   Colors are dynamic by position: violins 1-2 gray, 3-6 alternating light/dark blue,
#   7+ a red gradient (matches real usage: 2 control constructs, 4 lamina-associated, the
#   rest speckle-associated -- see project memory's intron-code decoder).
#
# This is a faithful rewrite of igoru's ~/R/trip/IntronTrip/plotSplicingEfficiencies.R.
# That script read `<name>.annotated_tables.txt` directly for exp.log10, which only works
# for the one legacy directory (IntronTRIP_combined2) where table_annotation.R had been run
# a second time per splicing-derived experiment, baking exp/exp.log10 straight into
# annotated_tables.txt -- a different, non-standard data layout from the rest of this
# pipeline (construct-level <mapName>.txt.annotated_tables.txt with NO exp columns, merged
# downstream against <expName>.combined_counts.txt by each analysis script, e.g. stage 12's
# annotateTRIP_positions_genuine_barcodes.R). Rewired here onto that same standard
# annotTable+countsTab merge (same exp.log10-preexisting guard and LP-style map_barcode
# handling as stage 12), so it works against IntronTRIP_combined/SON_combined-shaped data.
# Also dropped: an unrelated, broken "GC by intron" ggplot block (referenced an undefined
# `dfShow`, would have errored on execution -- copy-paste leftover from a different script,
# not part of this one) and an `n_samples <- nrow(dfShow)` line that would have crashed the
# same way before ever reaching the real violin plot; the dynamic (rather than hardcoded-16)
# sample-count color logic the original's own comment already called for is restored using
# the actual number of violins instead.
# =============================================================================
library(ggplot2)

path = "/home/projects/ulitsky/shared/data/human/TRIP/IntronTRIP_combined/"
expListFile = "splicing_efficiencies.exp_list.txt"
plot_name = "splicingEfficiencies.pdf"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs) > 0)
{
  path = myArgs[1]
  expListFile = myArgs[2]
  plot_name = myArgs[3]
}
print(path)
print(expListFile)
print(plot_name)

expList = read.table(paste0(path, expListFile), sep = "\t", header = FALSE,
                      col.names = c("expName", "niceName", "nameMap"),
                      colClasses = "character")

annotTables = list()
df_violin = list()
for (ind in seq_len(nrow(expList)))
{
  expName = expList$expName[ind]
  niceName = expList$niceName[ind]
  nameMap = expList$nameMap[ind]
  print(paste(expName, niceName, nameMap))

  if (is.null(annotTables[[nameMap]])) {
    annotTables[[nameMap]] = read.table(paste0(path, nameMap, ".annotated_tables.txt"), sep = "\t", quote = "")
  }
  annotTable = annotTables[[nameMap]]

  countsTab = read.table(paste0(path, expName, ".combined_counts.txt"), sep = "\t", quote = "")
  # Some combined_counts.txt files (e.g. *.spliced_unspliced.combined_counts.txt) already
  # carry their own exp.log10 -- a splicing ratio bounded [0,1], NA when undefined -- which
  # is a completely different metric from log10(1+exp). Track whether it was already
  # present so it doesn't get silently overwritten with the wrong formula (same guard as
  # stage 12/13/15/16's equivalent merges).
  hasOwnExpLog10 = "exp.log10" %in% names(countsTab)
  countsTab$Barcode = rownames(countsTab)
  annotTable$Barcode = rownames(annotTable)
  if ("map_barcode" %in% names(countsTab)) {
    # LP-style two-barcode data -- see stage 12's own comments for why this join is
    # keyed on map_barcode rather than a direct Barcode-to-Barcode merge.
    countsTab$recomb_barcode = countsTab$Barcode
    countsTab$Barcode = NULL
    merged = merge(annotTable, countsTab, by.x = c("Barcode","Chr","Pos","Strand"), by.y = c("map_barcode","Chr","Pos","Strand"))
  } else {
    merged = merge(annotTable, countsTab, by = c("Barcode","Chr","Pos","Strand"))
  }
  if (!hasOwnExpLog10) {
    merged$exp.log10 = log(1 + merged$exp, 10)
  }

  df_violin[[ind]] = data.frame(sample = niceName, exp.log10 = merged$exp.log10)
}
df_violin = do.call(rbind, df_violin)

# expListFile row order (first appearance of each niceName) is the plot order.
niceNamesOrdered = unique(expList$niceName)
df_violin$sample = factor(df_violin$sample, levels = niceNamesOrdered)

# Restrict to [0,1] before plotting -- exp.log10 here is a splicing-ratio fraction, not
# log10(1+exp).
df_violin = df_violin[
  !is.na(df_violin$exp.log10) &
    df_violin$exp.log10 >= 0 &
    df_violin$exp.log10 <= 1,
]

# Dynamic colors for however many violins are actually shown.
n_samples = length(niceNamesOrdered)
sample_colors = rep(NA, n_samples)
if (n_samples >= 1) sample_colors[1:min(2, n_samples)] = "gray70"
if (n_samples >= 3) {
  sample_colors[3:min(6, n_samples)] = rep(c("#9ecae1", "#08519c"), length.out = length(3:min(6, n_samples)))
}
if (n_samples >= 7) {
  sample_colors[7:n_samples] = colorRampPalette(c("#fcbba1", "#99000d"))(length(7:n_samples))
}
names(sample_colors) = niceNamesOrdered

pdf(plot_name, width = 12, useDingbats = F)
ggplot(df_violin, aes(x = sample, y = exp.log10, fill = sample)) +
  geom_violin(draw_quantiles = c(.25,.5,.75), quantile.linewidth = c(.3,.6,.3),
    trim = TRUE,
    bounds = c(0, 1),
    color = "black",
    scale = "width"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Splicing efficiency"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
