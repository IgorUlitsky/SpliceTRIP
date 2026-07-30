# =============================================================================
# Stage 12 — per-construct QC report: barcode overlap, ChromHMM/SPIN vs expression
# =============================================================================
# Usage
#   Rscript 12_annotateTRIP_positions_genuine_barcodes.R <path> <mapName>.txt <namesExp_file> <outFile.pdf>
#
# Inputs   <mapName>.txt.annotated_tables.txt   (stage 10)
#          <expName>.combined_counts.txt        (stage 07, or any pipeline producing the
#                                                same columns -- e.g. a splicing-ratio
#                                                exp.log10, auto-detected and preserved;
#                                                see examples/exp_lists/ for real ones)
#   namesExp_file: one <expName> per line (NOT the <mapName><TAB><expName> pairs stage 11
#                  uses -- e.g. examples/exp_lists/single_combined.exp_list.txt's 2nd
#                  column, filtered to one construct). run_pipeline.sh's "annotated"
#                  target builds this per construct from ### COMBOS automatically.
#
# Outputs  <outFile.pdf>                 barcode-overlap heatmap + per-library QC pages
#          <outFile.pdf>.summary.pdf     cross-library ChromHMM/SPIN summary heatmaps
#          <outFile.pdf>.summary.txt / .counts.txt / .summary_SPIN.txt / .counts_SPIN.txt
#          <outFile.pdf>.<expName>.raw_ChromHMM.txt / .raw_SPIN.txt   per-barcode
#                                        (state, value) pairs behind each boxplot, for
#                                        pooling across an intron's replicates elsewhere
#
# Run once PER construct (loop over every distinct ### FASTQ_MAPS mapName), unlike stage
# 11 which takes the whole construct list in one call.
#
# This is igoru's `annotateTRIP_positions_genuine_barcodes.R` verbatim, aside from an
# added doc header and dropping one line (`chromSizes <- read.table(".../hg38.chrom.
# sizes", ...)`) that loaded a value never actually referenced anywhere else in the
# script -- dead code, not a behaviour change. No other hardcoded external paths: unlike
# stages 08-10, everything else it reads is other pipeline stages' own output.
# =============================================================================
library(rtracklayer)
library(pheatmap)
library(ggplot2)
library(dplyr)

# Draws a ggplot boxplot from a named list of numeric vectors (same shape as the
# base-R `boxData` lists used throughout this script), one box per list element.
ggBoxplot <- function(boxData, labels = NULL, title = NULL, xlab = NULL, ylab = NULL, log_y = FALSE, angle = 0) {
  if (is.null(names(boxData))) names(boxData) <- seq_along(boxData)
  plot_df <- do.call(rbind, lapply(names(boxData), function(k) {
    data.frame(group = k, value = boxData[[k]])
  }))
  plot_df$group <- factor(plot_df$group, levels = names(boxData))
  if (is.null(labels)) labels <- paste0(names(boxData), "\n", sapply(boxData, length))
  p <- ggplot(plot_df, aes(x = group, y = value)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.shape = 20, outlier.size = 0.8) +
    scale_x_discrete(labels = labels) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw()
  if (log_y) p <- p + scale_y_log10()
  if (angle > 0) p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 0.5))
  p
}

# Blank divider page with just a title, used to mark the start of each library's section.
ggSectionHeader <- function(title) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title, size = 10, fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

mainChroms <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# Draws a ggplot pie chart from a named/table-like vector of counts (e.g. table(df$Chr)).
# If level_order is given, slices/legend follow that order instead of being sorted by count.
ggPie <- function(counts, title = NULL, level_order = NULL) {
  plot_df <- data.frame(category = names(counts), count = as.numeric(counts))
  if (is.null(level_order)) {
    plot_df$category <- factor(plot_df$category, levels = plot_df$category[order(plot_df$count)])
  } else {
    plot_df$category <- factor(plot_df$category, levels = level_order)
  }
  ggplot(plot_df, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.2) +
    coord_polar("y", start = 0) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

# Boxplot of expression by a categorical genomic state (ChromHMM or SPIN), ordered by
# mean, restricted to states with enough support (>=10 barcodes or >=2% of the library).
# Returns per-state mean/count, used to populate the cross-library summary matrices.
# Also saves the raw per-barcode (state, value) pairs behind the boxplot -- restricted to
# the same >=10/2%-supported states shown in the plot -- to <plot_name>.<name>.raw_<stateCol>.txt,
# so downstream scripts (analyzeIntronTRIPSPIN.R / analyzeIntronTRIPHMM.R) can pool these
# across an intron's replicate libraries and run their own Kruskal-Wallis test on the real
# distributions, rather than just on the already-collapsed per-state means saved elsewhere.
plotStateExpressionBoxplot <- function(df, stateCol, title, isP, plot_name, name) {
  yCol = if (isP) "exp.norm" else "exp.log10"
  yLab = if (isP) "Normalized Expression" else "Expression, log10(x+1) scale"
  sp = split(df, df[[stateCol]])
  sp = sp[sapply(sp, nrow) >= max(10, 0.02 * nrow(df))]
  myOrder = order(sapply(lapply(sp, "[[", yCol), mean))
  boxData = lapply(sp[myOrder], "[[", yCol)
  print(ggBoxplot(boxData, labels = paste0(names(boxData), " (", sapply(boxData, length), ")"),
                   title = title, ylab = yLab, angle = 90))
  rawDf = do.call(rbind, lapply(names(boxData), function(s) data.frame(state = s, value = boxData[[s]])))
  write.table(rawDf, paste0(plot_name, ".", name, ".raw_", stateCol, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(means = sapply(boxData, mean, na.rm = TRUE), counts = sapply(boxData, length))
}

# Two stacked bar charts crossing a categorical genomic state (ChromHMM or SPIN) against
# expression-level bins (df$exp.level, already computed -- shared between whichever state
# columns are being cross-tabulated, since it depends only on expression): state
# composition per expression level, and expression-level composition per state.
plotStateVsExpressionLevel <- function(df, stateCol, statesLabel, name, isP, expLevelColors) {
  levelLab = if (isP) "Normalized Expression Level" else "Log10 (Expression + 1) Level"
  titlePrefix = if (isP) "Normalized Expression Levels by" else "Expression Levels by"

  byState <- df %>%
    filter(!is.na(.data[[stateCol]])) %>%
    group_by(exp.level, across(all_of(stateCol))) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(across(all_of(stateCol))) %>%
    mutate(percentage = count / sum(count) * 100)

  p1 <- ggplot(byState, aes(x = .data[[stateCol]], y = percentage, fill = exp.level)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = "black", size = 2) +
    labs(title = paste(titlePrefix, statesLabel, name), x = statesLabel, y = "Percentage", fill = levelLab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = expLevelColors)
  print(p1)

  byLevel <- df %>%
    filter(!is.na(.data[[stateCol]])) %>%
    group_by(across(all_of(stateCol)), exp.level) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(exp.level) %>%
    mutate(percentage = count / sum(count) * 100)

  custom_colors <- c("#f44336", "#e81e63", "#9c27b0", "#673ab7", "#3f51b5",
                     "#2196f3", "#03a9f4", "#00bcd4", "#009688", "#4caf50",
                     "#8bc34a", "#cddc39", "#ffeb3b", "#ffc107", "#ff9800",
                     "#ff5722", "#ffffff", "#000000")

  p2 <- ggplot(byLevel, aes(x = exp.level, y = percentage, fill = .data[[stateCol]])) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = "black", size = 2) +
    labs(title = paste(statesLabel, "Distribution by Expression Level", name), x = "Expression Level", y = "Percentage", fill = statesLabel) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = custom_colors)
  print(p2)
}

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined//"
nameMap <- "5ss.txt"
namesExp <- read.table("/home/projects/ulitsky/shared/data/human/TRIP/single_combined/names.5ss.txt")$V1

# Arguments from the terminal
myArgs = commandArgs(trailingOnly = T)
plot_name="temp.pdf"
if (length(myArgs)>0)
{
  path = myArgs[1]
  nameMap = myArgs[2]
  namesExp_file = myArgs[3]
  plot_name = myArgs[4]
  namesExp<-read.table(paste0(path, namesExp_file))$V1
}
print(path)
print(namesExp)
print(plot_name)


niceNames = gsub("_Expression", "", namesExp)
niceNames = gsub("Exp_", "", niceNames)
niceNames = gsub("_S[0-9].*", "", niceNames)

annotTable =  read.table(paste0(path, nameMap, ".annotated_tables.txt"), sep = "\t", quote = "")

allBarcodes = c()
dfs = list()
for (ind in 1:length(namesExp)) {
  print(ind)
  nameExp = namesExp[ind]
  name = niceNames[ind]
  countsTab =  read.table(paste0(path, nameExp, ".combined_counts.txt"), sep = "\t", quote = "")
  # Some combined_counts.txt files (e.g. *.spliced_unspliced.combined_counts.txt) already
  # carry their own exp.log10 -- a splicing ratio bounded [0,1], NA when undefined (e.g.
  # exp_total==0) -- which is a completely different metric from log10(1+exp) below. Track
  # whether it was already present so it doesn't get silently overwritten with the wrong
  # formula for those files (confirmed on real data: this produced e.g. exp.log10=1.10 for
  # a metric that should never exceed 1).
  hasOwnExpLog10 = "exp.log10" %in% names(countsTab)
  countsTab$Barcode = rownames(countsTab)
  annotTable$Barcode = rownames(annotTable)
  if ("map_barcode" %in% names(countsTab)) {
    # LP-style two-barcode data: countsTab's own row identity is the downstream
    # recomb_barcode -- an independent reporter measurement per recombination event,
    # multiple of which can share one genomic integration -- living in a completely
    # different barcode space (18nt) than annotTable's map-barcode-keyed Barcode
    # (16nt), so a direct Barcode-to-Barcode merge always returns 0 rows. Join on
    # map_barcode instead, but keep recomb_barcode as the row identity throughout
    # (matching this script's "one row = one physical barcode" semantics, and the
    # recomb_barcode/map_barcode naming already used in normalizeTRIP_LP.R).
    countsTab$recomb_barcode = countsTab$Barcode
    countsTab$Barcode = NULL
    merged = merge(annotTable, countsTab, by.x = c("Barcode","Chr","Pos","Strand"), by.y = c("map_barcode","Chr","Pos","Strand"))
    rownames(merged) = merged$recomb_barcode
  } else {
    merged = merge(annotTable,countsTab,by=c("Barcode","Chr","Pos","Strand"))
    rownames(merged)=merged$Barcode
  }
  if (!hasOwnExpLog10) {
    merged$exp.log10 = log(1+merged$exp,10)
  }
  merged$exp.norm = log((1 + merged$exp) / (1 + merged$norm), 2)
  merged$Row.names=NULL
  dfs[[ind]]=merged
  dfs[[ind]]$CompositePos = paste0(dfs[[ind]]$Chr, ":", dfs[[ind]]$Pos)
  allBarcodes = c(allBarcodes, rownames(dfs[[ind]]))
}
tBarcodes = table(allBarcodes)
sort(tBarcodes, decreasing = T)[1:20]

truncate_row_names <- function(df) {
  rownames(df) <- substring(rownames(df), 1, 16)
  return(df)
}

#dfs_trunc <- lapply(dfs, truncate_row_names)
dfs_trunc <-dfs
oMat = matrix(nrow = length(dfs_trunc), ncol = length(dfs_trunc), NA)
for (i in 1:length(dfs_trunc)) {
  for (j in 1:length(dfs_trunc)) {
    oMat[i, j] = length(intersect(rownames(dfs_trunc[[i]]), rownames(dfs_trunc[[j]])))
  }
}

pdf(plot_name, width = 12)

dfO = data.frame(oMat)
names(dfO) = niceNames
rownames(dfO) = niceNames

pheatmap(log(1 + dfO, 10),main="Overlap between barcodes (log10)")


# Vector with datasets that should be normalized (with promoter)
names_with_P_list = c('5ss_Expression_rep1_S56', '5ss_Expression_rep2_S57',
                      'CPA_Expression_rep2_S54', 'CPA_Expression_rep3_S55')

allStates = c()
for (ind in 1:length(namesExp)) {
  df = dfs[[ind]]
  allStates = unique(c(allStates,na.omit(df$ChromHMM)))
}

allSPINStates = c()
for (ind in 1:length(namesExp)) {
  df = dfs[[ind]]
  allSPINStates = unique(c(allSPINStates,na.omit(df$SPIN)))
}

dfAllChromHMM = data.frame(row.names=allStates)
dfAllCounts = data.frame(row.names=allStates)
dfAllSPIN = data.frame(row.names=allSPINStates)
dfAllSPINCounts = data.frame(row.names=allSPINStates)
for (ind in 1:length(namesExp)) {
  print(ind)
  nameExp = namesExp[ind]
  name = nameExp
  print(paste(ind, name))
  df = dfs[[ind]]
  print(ggSectionHeader(name))

  tPos = table(paste0(df$CompositePos))

  boxData = list(df[ (df$CompositePos %in% names(tPos)[tPos > 1]), "norm"] + 1,
                 df[df$CompositePos %in% names(tPos)[tPos == 1], "norm"] + 1)
  names(boxData) = c("Multi\nintegrants", "Single\nintegrants")
  print(ggBoxplot(boxData, title = name, ylab = "1+(norm reads)", log_y = TRUE))

  sp = split(df, tBarcodes[rownames(df)])
  boxData = lapply(lapply(sp, "[[", "norm"), "+", 1)
  print(ggBoxplot(boxData, title = name,
                   xlab = "Number of times barcode appears in different libraries",
                   ylab = "norm", log_y = TRUE))

  if (sum(grepl("chr11:522674", df$CompositePos)) > 0) {
    par(mfrow = c(1, 2))
    heatscatter(1 + df[grepl("chr11:522674", df$CompositePos), "norm"],
                1 + df[grepl("chr11:522674", df$CompositePos), "exp"], log = "xy",
                main = "Globin", xlab = "Normalizer", ylab = "Expression")
    heatscatter(1 + df[!grepl("chr11:522674", df$CompositePos), "norm"],
                1 + df[!grepl("chr11:522674", df$CompositePos), "exp"], log = "xy",
                main = "Non-globin", xlab = "Normalizer", ylab = "Expression")
  }

  #NORM FILTER
  df = subset(df, df$norm >= 10)

  names(df)
  chrCounts = table(df$Chr)
  chrCounts = chrCounts[names(chrCounts) %in% mainChroms]
  print(ggPie(chrCounts, title = name, level_order = mainChroms[mainChroms %in% names(chrCounts)]))

  isP = name %in% names_with_P_list

  # ChromHMM states vs expression
  resChromHMM = plotStateExpressionBoxplot(df, "ChromHMM", name, isP, plot_name, name)
  dfAllChromHMM[,name]=resChromHMM$means[allStates]
  dfAllCounts[,name]=resChromHMM$counts[allStates]
  dfAllCounts[,name][is.na(dfAllCounts[,name])]=0

  # SPIN states vs expression
  resSPIN = plotStateExpressionBoxplot(df, "SPIN", paste(name, "- SPIN"), isP, plot_name, name)
  dfAllSPIN[,name]=resSPIN$means[allSPINStates]
  dfAllSPINCounts[,name]=resSPIN$counts[allSPINStates]
  dfAllSPINCounts[,name][is.na(dfAllSPINCounts[,name])]=0

  # Chr vs expression
  if (name %in% names_with_P_list) {
    sp = split(df, df$Chr)
    sp = sp[sapply(sp, nrow) >= 10]
    myOrder = order(sapply(lapply(sp, "[[", "exp.norm"), mean))
    boxData = lapply(sp[myOrder], "[[", "exp.norm")
    print(ggBoxplot(boxData, labels = paste0(names(boxData), " (", sapply(boxData, length), ")"),
                     title = name, ylab = "Normalized Expression", angle = 90))
  } else {
    sp = split(df, df$Chr)
    sp = sp[sapply(sp, nrow) >= 10]
    myOrder = order(sapply(lapply(sp, "[[", "exp.log10"), mean))
    boxData = lapply(sp[myOrder], "[[", "exp.log10")
    print(ggBoxplot(boxData, labels = paste0(names(boxData), " (", sapply(boxData, length), ")"),
                     title = name, ylab = "Expression, log10(x+1) scale", angle = 90))
  }

  # Expression-level (quantile) bins -- shared between the ChromHMM and SPIN
  # cross-tabulations below, computed once per library since it depends only on
  # expression, not on either state annotation.

  create_unique_breaks <- function(values, num_quantiles) {
    repeat {
      quantile_breaks <- quantile(values, probs = seq(0, 1, length.out = num_quantiles + 1), na.rm = TRUE)
      breaks <- unique(c(-Inf, 0, quantile_breaks, Inf))
      if (length(breaks) == length(unique(breaks))) {
        return(breaks)
      }
      num_quantiles <- num_quantiles - 1
      if (num_quantiles < 1) {
        stop("Unable to create unique breaks with less than 2 quantiles.")
      }
    }
  }

  if (name %in% names_with_P_list) {
    breaks <- c(-Inf, 0, quantile(df$exp.norm[df$exp.norm > 0], seq(0.1, 1, by = 0.1)), Inf)
    labels <- c("Below 0 norm.exp", paste0("Quantile ", 1:11,
                                           " (", round(breaks[-c(1, length(breaks), length(breaks) - 1)], 2),
                                           " - ", round(breaks[-c(1:2, length(breaks))], 2), ")"))

    df$exp.level <- cut(df$exp.norm,
                        breaks = breaks,
                        labels = labels,
                        include.lowest = TRUE)

    expLevelColors <- c("Below 0 norm.exp" = "grey",
                         setNames(scales::seq_gradient_pal("white", "blue")(seq(0, 1, length.out = 11)), labels[-1]))
  } else{
    num_quantiles <- 6

    breaks <- create_unique_breaks(df$exp.log10[df$exp.log10 > 0], num_quantiles)

    labels <- c("Zero Exp",
                paste0("Quantile ", seq_len(length(breaks) - 2),
                       " (", round(breaks[2:(length(breaks) - 2)], 2),
                       " - ", round(breaks[3:(length(breaks) - 1)], 2), ")"))

    df$exp.level <- cut(df$exp.log10,
                        breaks = breaks,
                        labels = labels,
                        include.lowest = TRUE)

    expLevelColors <- setNames(scales::seq_gradient_pal("white", "blue")(seq(0, 1, length.out = 10)), labels)
  }

  # ChromHMM state vs Expression level (both directions)
  plotStateVsExpressionLevel(df, "ChromHMM", "ChromHMM", name, isP, expLevelColors)

  # SPIN state vs Expression level (both directions)
  plotStateVsExpressionLevel(df, "SPIN", "SPIN", name, isP, expLevelColors)
}

dev.off()

pdf(paste0(plot_name,".summary.pdf"),width = 12)
write.table(dfAllChromHMM,paste0(plot_name,".summary.txt"),sep="\t",quote=F)
write.table(dfAllCounts,paste0(plot_name,".counts.txt"),sep="\t",quote=F)
pheatmap(na.omit(dfAllChromHMM))
pheatmap(dfAllChromHMM,cluster_rows = F,cluster_cols = F)

write.table(dfAllSPIN,paste0(plot_name,".summary_SPIN.txt"),sep="\t",quote=F)
write.table(dfAllSPINCounts,paste0(plot_name,".counts_SPIN.txt"),sep="\t",quote=F)
pheatmap(na.omit(dfAllSPIN))
pheatmap(dfAllSPIN,cluster_rows = F,cluster_cols = F)
dev.off()
