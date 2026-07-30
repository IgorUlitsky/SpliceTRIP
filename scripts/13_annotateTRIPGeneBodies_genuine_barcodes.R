# =============================================================================
# Stage 13 — per-construct gene-body positional analysis (GC/speckle-distance/host-
# expression/TSS-TES-intron-distance correlations, density-terrain scatters)
# =============================================================================
# Usage
#   Rscript 13_annotateTRIPGeneBodies_genuine_barcodes.R <path> <mapName>.txt <namesExp_file> <outFile.pdf> [speckleDir]
#
# Inputs   <mapName>.txt.annotated_tables.txt   (stage 10)
#          <expName>.combined_counts.txt        (stage 07, or any pipeline producing the
#                                                same columns -- e.g. a splicing-ratio
#                                                exp.log10, auto-detected and preserved;
#                                                see examples/exp_lists/ for real ones)
#   namesExp_file: one <expName> per line, same format/derivation as stage 12's.
#
# Outputs  <outFile.pdf>                per-library correlation scatters (heatscatter +
#                                        ggplot density-terrain) against GC content,
#                                        speckle distance, host gene expression,
#                                        TSS/TES/intron-boundary distance, relative
#                                        position, intron length
#          <outFile.pdf>summary.pdf     cross-library correlation summary (bar charts +
#                                        heatmaps)
#          <outFile.pdf>.summary.txt / .counts.txt / .pvalue.xt   the underlying Spearman
#                                        rho/n/p values behind every correlation above
#
# Run once PER construct (loop over every distinct ### FASTQ_MAPS mapName), same as
# stage 12 -- unlike stage 11 which takes the whole construct list in one call.
#
# This is igoru's `annotateTRIPGeneBodies_genuine_barcodes.R` verbatim, aside from an
# added doc header and its one hardcoded external path (the MANE-to-speckle-distance
# reference directory) turned into an optional 5th CLI arg (keeps the original hardcoded
# path as default when omitted, same length-guarded pattern as stage 06's barcodeLen/
# dpniiBedPath) so the config drives it instead.
# =============================================================================
library(LSD)
library(pheatmap)
library(ggplot2)
library(ggpubr)

# ggplot filled-contour ("terrain-style") 2D density plot, shown alongside each heatscatter().
ggDensityTerrain <- function(x, y, xlab = NULL, ylab = NULL, title = NULL, log_x = FALSE) {
  df <- data.frame(x = x, y = y)
  df <- df[is.finite(df$x) & is.finite(df$y), ]
  # Match heatscatter(..., log="x")'s convention: it can't log-transform x<=0, so it
  # drops those rows before computing its displayed correlation. Do the same here so
  # (a) our R/P matches heatscatter's for the same data, and (b) our own R/P matches
  # what's actually drawn -- scale_x_log10() below would silently drop these points
  # from the plot anyway, so leaving them in the correlation was already inconsistent
  # with our own plot.
  if (log_x) df <- df[df$x > 0, ]
  ct <- tryCatch(cor.test(df$x, df$y, method = "spearman"), error = function(e) NULL)
  corLabel <- if (!is.null(ct)) sprintf("R = %3.2g, P = %3.2g", ct$estimate, ct$p.value) else "R = NA, P = NA"
  # MASS::kde2d (used by stat_density_2d) needs a strictly positive bandwidth in both
  # dimensions, or it errors out with "bandwidths must be strictly positive". A column
  # can look fine overall yet collapse to zero IQR within one particular x/y subset
  # (e.g. exp.log10 is >=50% zero within some subsets), so pre-checking unique-value
  # counts isn't reliable -- check the actual bandwidth directly instead. This is the
  # same check MASS::kde2d does internally, so it's cheap (no wasted KDE fit) and exact.
  canFitDensity <- nrow(df) >= 3 &&
    isTRUE(tryCatch(MASS::bandwidth.nrd(df$x) > 0, error = function(e) FALSE)) &&
    isTRUE(tryCatch(MASS::bandwidth.nrd(df$y) > 0, error = function(e) FALSE))
  p <- ggplot(df, aes(x = x, y = y))
  if (canFitDensity) {
    p <- p +
      stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, bins = 25, n = 50) +
      scale_fill_distiller(palette = "Greens", direction = 1, name = "density")
  }
  # Cap the point overlay so PDF vector rendering doesn't choke on very dense libraries;
  # density/correlation above are still computed on the full data.
  pointDf <- if (nrow(df) > 4000) df[sample.int(nrow(df), 4000), ] else df
  p <- p +
    geom_point(data = pointDf, color = "black", size = 0.6, alpha = 0.4) +
    labs(title = title, subtitle = corLabel, x = xlab, y = ylab) +
    theme_bw(base_size = 16)
  if (log_x) p <- p + scale_x_log10()
  p
}

# Blank divider page with just a title, used to mark the start of each library's section.
ggSectionHeader <- function(title) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title, size = 10, fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined/"

namesExp<-read.table("/home/projects/ulitsky/shared/data/human/TRIP/single_combined/names.5ss.txt")$V1

plot_name <- '/home/projects/ulitsky/shared/data/human/TRIP/single_combined/annotateTRIPGeneBodies_5ss.pdf'
speckleDir <- "/home/projects/ulitsky/shared/data/human/speckles/"
# Arguments from the terminal
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  nameMap = myArgs[2]
  namesExp_file = myArgs[3]
  plot_name = myArgs[4]
  # speckleDir (5th arg) is optional -- keeps its default above when omitted, same
  # length-guarded pattern as stage 06's barcodeLen/dpniiBedPath.
  if (length(myArgs)>=5 && nzchar(myArgs[5])) speckleDir = myArgs[5]
  namesExp<-read.table(paste0(path, namesExp_file))$V1
}
print(path)
print(namesExp)
print(plot_name)

niceNames = gsub("_Expression", "", namesExp)
niceNames = gsub("Exp_", "", niceNames)
niceNames = gsub("_S[0-9].*", "", niceNames)
niceNames = gsub("rep([0-9]+)", "\\1", niceNames, ignore.case = TRUE)
niceNames = gsub("pool([0-9]+)", "\\1", niceNames, ignore.case = TRUE)



#names_with_P_list = c('5ss_Expression_rep1_S56', '5ss_Expression_rep2_S57',
#                      'CPA_Expression_rep2_S54', 'CPA_Expression_rep3_S55')
names_with_P_list = c()

speckleNames = read.table(paste0(speckleDir, "speckle_names.txt"))$V1
maneSpeckleData=list()
for (speckleName in speckleNames)
{
  maneSpeckleData[[speckleName]]=read.table(gzfile(paste0(speckleDir, "MANE.",speckleName,".tab.gz")),header=F,sep="\t",row.names=1)
}

avSpeckleDistance = data.frame(row.names=rownames(maneSpeckleData[[1]]),avSpeckleDistance=apply(do.call("cbind",lapply(maneSpeckleData,"[[","V5"))[,c("GSE81553_Speckle_distance_Condition1","GSE81553_Speckle_distance_Condition2","GSE81553_Speckle_distance_Condition3")],1,mean,na.rm=T))

annotTable =  read.table(paste0(path, nameMap, ".annotated_tables.txt"), sep = "\t", quote = "")


pdf(plot_name, useDingbats = F, width = 14)
tabs = list()
dfSum = data.frame(row.names=namesExp)
dfP = data.frame(row.names=namesExp)
dfCounts = data.frame(row.names=namesExp)
for (name in namesExp)
{
  print(name)
  countsTab =  read.table(paste0(path, name, ".combined_counts.txt"), sep = "\t", quote = "")
  # Some combined_counts.txt files (e.g. *.spliced_unspliced.combined_counts.txt) already
  # carry their own exp.log10 -- a splicing ratio bounded [0,1], NA when undefined -- which
  # is a completely different metric from log10(1+exp) below. Track whether it was already
  # present so it doesn't get silently overwritten with the wrong formula for those files
  # (confirmed elsewhere in this pipeline: an identical unconditional overwrite in
  # annotateTRIP_positions_genuine_barcodes.R produced values >1 for a metric that should
  # never exceed 1).
  hasOwnExpLog10 = "exp.log10" %in% names(countsTab)
  countsTab$Barcode = rownames(countsTab)
  annotTable$Barcode = rownames(annotTable)
  merged = merge(annotTable,countsTab,by=c("Barcode","Chr","Pos","Strand"))
  if (!hasOwnExpLog10) {
    merged$exp.log10 = log(1+merged$exp,10)
  }
  merged$exp.norm = log((1 + merged$exp) / (1 + merged$norm), 2)
  rownames(merged)=merged$Barcode
  merged$Row.names=NULL

  tabs[[name]] = merged
  print(ggSectionHeader(name))

  tabs[[name]][["avSpeckleDistance"]]=avSpeckleDistance[tabs[[name]]$MANE_all,"avSpeckleDistance"]

  # Boxplot: Normalized expression / log10(exp + 1)
  counts <- c(
    sum(!is.na(tabs[[name]]$MANE_all)),
    sum(is.na(tabs[[name]]$MANE_all)),
    sum(!is.na(tabs[[name]]$MANE_exonic)),
    sum(is.na(tabs[[name]]$MANE_exonic)),
    sum(!is.na(tabs[[name]]$MANE_antisense)),
    sum(is.na(tabs[[name]]$MANE_antisense)),
    sum(!is.na(tabs[[name]]$MANE_exonic_antisense)),
    sum(is.na(tabs[[name]]$MANE_exonic_antisense)),
    sum(!is.na(tabs[[name]]$gencode_noncoding))
  )

  labels_boxplot <- paste(c("Sense\nall", "Not Sense\nall",
                            "Sense\nexonic", "Not Sense\nexonic",
                            "Antisense\nall", "Not Antisense\nall",
                            "Antisense\nexonic", "Not Antisense\nexonic",
                            "Gencode sense\nNoncoding"),
                          "(", counts, ")", sep = "")
  valueCol <- if (name %in% names_with_P_list) "exp.norm" else "exp.log10"
  ylabText <- if (name %in% names_with_P_list) "Normalized expression" else "Log10 (Expression + 1)"

  boxGroups <- list(
    tabs[[name]][!is.na(tabs[[name]]$MANE_all), valueCol],
    tabs[[name]][is.na(tabs[[name]]$MANE_all), valueCol],
    tabs[[name]][!is.na(tabs[[name]]$MANE_exonic), valueCol],
    tabs[[name]][is.na(tabs[[name]]$MANE_exonic), valueCol],
    tabs[[name]][!is.na(tabs[[name]]$MANE_antisense), valueCol],
    tabs[[name]][is.na(tabs[[name]]$MANE_antisense), valueCol],
    tabs[[name]][!is.na(tabs[[name]]$MANE_exonic_antisense), valueCol],
    tabs[[name]][is.na(tabs[[name]]$MANE_exonic_antisense), valueCol],
    tabs[[name]][!is.na(tabs[[name]]$gencode_noncoding), valueCol]
  )
  # "Sense"/"Not Sense"/"Gencode sense" groups vs "Antisense"/"Not Antisense" groups
  strandType <- c("Sense", "Sense", "Sense", "Sense", "Antisense", "Antisense", "Antisense", "Antisense", "Sense")

  boxplot_df <- do.call(rbind, lapply(seq_along(boxGroups), function(i) {
    data.frame(group = labels_boxplot[i], value = boxGroups[[i]], strand = strandType[i])
  }))
  boxplot_df$group <- factor(boxplot_df$group, levels = labels_boxplot)

  p <- ggplot(boxplot_df, aes(x = group, y = value, fill = strand)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 20, outlier.size = 0.8) +
    scale_fill_manual(values = c(Sense = "steelblue", Antisense = "indianred"), name = NULL) +
    labs(title = name, x = NULL, y = ylabText) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  print(p)
  # Heatscatter: MANE_all_exp.norm vs exp.norm / log10(exp + 1)
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$GC, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Local GC content)",
                ylab = "Integration normalized expression",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations"))
    print(ggDensityTerrain(tabs[[name]]$GC, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Local GC content", ylab = "Integration normalized expression",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations")))
  } else{
    dfSum[name,"CorrGC"]=cor.test(tabs[[name]]$GC, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrGC"]=cor.test(tabs[[name]]$GC, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrGC"]=sum(!is.na(tabs[[name]]$GC) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$GC, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Local GC content",
                ylab = "Integration log10(expression + 1)",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations"))
    print(ggDensityTerrain(tabs[[name]]$GC, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Local GC content", ylab = "Integration log10(expression + 1)",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations")))
  }

  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$avSpeckleDistance, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Gene distance to speckles)",
                ylab = "Integration normalized expression",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations"))
    print(ggDensityTerrain(tabs[[name]]$avSpeckleDistance, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Gene distance to speckles", ylab = "Integration normalized expression",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations")))
  } else{
    dfSum[name,"CorrSpeckleDistance"]=cor.test(tabs[[name]]$avSpeckleDistance, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrSpeckleDistance"]=cor.test(tabs[[name]]$avSpeckleDistance, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrSpeckleDistance"]=sum(!is.na(tabs[[name]]$avSpeckleDistance) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$avSpeckleDistance, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Gene distance to speckles",
                ylab = "Integration log10(expression + 1)",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations"))
    print(ggDensityTerrain(tabs[[name]]$avSpeckleDistance, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Gene distance to speckles", ylab = "Integration log10(expression + 1)",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations")))
  }

  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$MANE_all_exp.norm, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Exonic Expression of the host gene (PolyAplus)",
                ylab = "Integration normalized expression",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations"))
    print(ggDensityTerrain(tabs[[name]]$MANE_all_exp.norm, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Exonic Expression of the host gene (PolyAplus)", ylab = "Integration normalized expression",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations")))
  } else{
    dfSum[name,"CorrExpHost"]=cor.test(tabs[[name]]$MANE_all_exp.norm, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpHost"]=cor.test(tabs[[name]]$MANE_all_exp.norm, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpHost"]=sum(!is.na(tabs[[name]]$MANE_all_exp.norm) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$MANE_all_exp.norm, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Exonic Expression of the host gene (PolyAplus)",
                ylab = "Integration log10(expression + 1)",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations"))
    print(ggDensityTerrain(tabs[[name]]$MANE_all_exp.norm, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Exonic Expression of the host gene (PolyAplus)", ylab = "Integration log10(expression + 1)",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_all)), "sense integrations")))
  }

  # Heatscatter: MANE_exonic_exp.norm vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$MANE_exonic_exp.norm, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Exonic Expression of the host gene (PolyAplus)",
                ylab = "Integration normalized expression",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_exonic)), "sense exonic integrations"))
    print(ggDensityTerrain(tabs[[name]]$MANE_exonic_exp.norm, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Exonic Expression of the host gene (PolyAplus)", ylab = "Integration normalized expression",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_exonic)), "sense exonic integrations")))
  }else{
    dfSum[name,"CorrExpHostExonic"]=cor.test(tabs[[name]]$MANE_exonic_exp.norm, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpHostExonic"]=cor.test(tabs[[name]]$MANE_exonic_exp.norm, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpHostExonic"]=sum(!is.na(tabs[[name]]$MANE_exonic_exp.norm) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$MANE_exonic_exp.norm, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Exonic Expression of the host gene (PolyAplus)",
                ylab = "Integration log10(expression + 1)",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_exonic)), "sense exonic integrations"))
    print(ggDensityTerrain(tabs[[name]]$MANE_exonic_exp.norm, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Exonic Expression of the host gene (PolyAplus)", ylab = "Integration log10(expression + 1)",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_exonic)), "sense exonic integrations")))
  }

  # Heatscatter: MANE_antisense_exp.norm vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$MANE_antisense_exp.norm, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Exonic Expression of the host gene (PolyAplus)",
                ylab = "Integration normalized expression",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_antisense)), "antisense integrations"))
    print(ggDensityTerrain(tabs[[name]]$MANE_antisense_exp.norm, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Exonic Expression of the host gene (PolyAplus)", ylab = "Integration normalized expression",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_antisense)), "antisense integrations")))
  }else{
    dfSum[name,"CorrExpHostExonicAntisense"]=cor.test(tabs[[name]]$MANE_antisense_exp.norm, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpHostExonicAntisense"]=cor.test(tabs[[name]]$MANE_antisense_exp.norm, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpHostExonicAntisense"]=sum(!is.na(tabs[[name]]$MANE_antisense_exp.norm) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$MANE_antisense_exp.norm, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Exonic Expression of the host gene (PolyAplus)",
                ylab = "Integration log10(expression + 1)",
                main = paste(name, sum(!is.na(tabs[[name]]$MANE_antisense)), "antisense integrations"))
    print(ggDensityTerrain(tabs[[name]]$MANE_antisense_exp.norm, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Exonic Expression of the host gene (PolyAplus)", ylab = "Integration log10(expression + 1)",
                title = paste(name, sum(!is.na(tabs[[name]]$MANE_antisense)), "antisense integrations")))
  }


  # Heatscatter: DistToTSS vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$DistToTSS, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Distance to TSS (nt)",
                ylab = "Normalized expression",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToTSS, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Distance to TSS (nt)", ylab = "Normalized expression", title = name))
  }else{
    dfSum[name,"CorrExpDistToTSS"]=cor.test(tabs[[name]]$DistToTSS, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpDistToTSS"]=cor.test(tabs[[name]]$DistToTSS, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpDistToTSS"]=sum(!is.na(tabs[[name]]$DistToTSS) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$DistToTSS, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Distance to TSS (nt)",
                ylab = "Log10(expression + 1)",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToTSS, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Distance to TSS (nt)", ylab = "Log10(expression + 1)", title = name))
  }

  # Heatscatter: DistToTES vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$DistToTES, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Distance to TES (nt)",
                ylab = "Normalized expression",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToTES, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Distance to TES (nt)", ylab = "Normalized expression", title = name))
  } else{
    dfSum[name,"CorrExpDistToTES"]=cor.test(tabs[[name]]$DistToTES, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpDistToTES"]=cor.test(tabs[[name]]$DistToTES, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpDistToTES"]=sum(!is.na(tabs[[name]]$DistToTES) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$DistToTES, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Distance to TES (nt)",
                ylab = "Log10(expression + 1)",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToTES, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Distance to TES (nt)", ylab = "Log10(expression + 1)", title = name))
  }

  # Heatscatter: RelPos vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$RelPos, tabs[[name]]$exp.norm, cor = T,
                xlab = "Relative position in the gene",
                ylab = "Normalized expression",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$RelPos, tabs[[name]]$exp.norm,
                xlab = "Relative position in the gene", ylab = "Normalized expression", title = name))
  } else{
    dfSum[name,"CorrExpRelPosGene"]=cor.test(tabs[[name]]$RelPos, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpRelPosGene"]=cor.test(tabs[[name]]$RelPos, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpRelPosGene"]=sum(!is.na(tabs[[name]]$RelPos) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$RelPos, tabs[[name]]$exp.log10, cor = T,
                xlab = "Relative position in the gene",
                ylab = "Log10(expression + 1)",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$RelPos, tabs[[name]]$exp.log10,
                xlab = "Relative position in the gene", ylab = "Log10(expression + 1)", title = name))
  }


  # Heatscatter: RelPos vs exp.norm for highly expressed hosts
  if (sum(!is.na(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$RelPos)) >= 10) {
    if (name %in% names_with_P_list) {
      heatscatter(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.norm, cor = T,
                  xlab = "Relative position in the gene",
                  ylab = "Normalized expression",
                  main = paste("Highly Expressed hosts", name))
      print(ggDensityTerrain(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.norm,
                  xlab = "Relative position in the gene", ylab = "Normalized expression",
                  title = paste("Highly Expressed hosts", name)))
    } else {
      dfSum[name,"CorrExpRelPosGeneHighlyExpressed"]=cor.test(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos,
                                                              tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10,method="spearman")$estimate
      dfP[name,"CorrExpRelPosGeneHighlyExpressed"]=cor.test(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos,
                                                              tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10,method="spearman")$p.value

      dfCounts[name,"CorrExpRelPosGeneHighlyExpressed"]=sum(!is.na(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos) & !is.na(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10),na.rm=T)

      heatscatter(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10, cor = T,
                  xlab = "Relative position in the gene",
                  ylab = "Log10(expression + 1)",
                  main = paste("Highly Expressed hosts", name))
      print(ggDensityTerrain(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$RelPos,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10,
                  xlab = "Relative position in the gene", ylab = "Log10(expression + 1)",
                  title = paste("Highly Expressed hosts", name)))
    }
  }

  # Heatscatter: DistToIntron5 vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$DistToIntron5, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Distance to intron 5' (nt)",
                ylab = "Normalized expression",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToIntron5, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Distance to intron 5' (nt)", ylab = "Normalized expression", title = name))
  } else {
    dfSum[name,"CorrExpDistToIntron5"]=cor.test(tabs[[name]]$DistToIntron5, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpDistToIntron5"]=cor.test(tabs[[name]]$DistToIntron5, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpDistToIntron5"]=sum(!is.na(tabs[[name]]$DistToIntron5) & !is.na(tabs[[name]]$exp.log10),na.rm=T)

    heatscatter(tabs[[name]]$DistToIntron5, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Distance to intron 5' (nt)",
                ylab = "Log10(expression + 1)",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToIntron5, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Distance to intron 5' (nt)", ylab = "Log10(expression + 1)", title = name))
  }

  # Heatscatter: DistToIntron3 vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$DistToIntron3, tabs[[name]]$exp.norm, log = "x", cor = T,
                xlab = "Distance to intron 3' (nt)",
                ylab = "Normalized expression",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToIntron3, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Distance to intron 3' (nt)", ylab = "Normalized expression", title = name))
  } else{
    dfSum[name,"CorrExpDistToIntron3"]=cor.test(tabs[[name]]$DistToIntron3, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpDistToIntron3"]=cor.test(tabs[[name]]$DistToIntron3, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpDistToIntron3"]=sum(!is.na(tabs[[name]]$DistToIntron3) & !is.na(tabs[[name]]$exp.log10),na.rm=T)

    heatscatter(tabs[[name]]$DistToIntron3, tabs[[name]]$exp.log10, log = "x", cor = T,
                xlab = "Distance to intron 3' (nt)",
                ylab = "Log10(expression + 1)",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$DistToIntron3, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Distance to intron 3' (nt)", ylab = "Log10(expression + 1)", title = name))
  }

  # Heatscatter: DistToIntron3 vs exp.norm for highly expressed hosts
  if (sum(!is.na(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length)) >= 10) {
    if (name %in% names_with_P_list) {
      heatscatter(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.norm, cor = T,
                  xlab = "Distance to intron 3'",
                  ylab = "Normalized expression",
                  main = paste("Highly Expressed hosts", name), log = "x")
      print(ggDensityTerrain(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.norm, log_x = TRUE,
                  xlab = "Distance to intron 3'", ylab = "Normalized expression",
                  title = paste("Highly Expressed hosts", name)))
    } else{
      dfSum[name,"CorrExpDistToIntron3HighHosts"]=cor.test(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3,
                                                           tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10,method="spearman")$estimate
      dfP[name,"CorrExpDistToIntron3HighHosts"]=cor.test(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3,
                                                           tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10,method="spearman")$p.value
      dfCounts[name,"CorrExpDistToIntron3HighHosts"]=sum(!is.na(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3) & !is.na(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10),na.rm=T)

      heatscatter(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10, cor = T,
                  xlab = "Distance to intron 3'",
                  ylab = "Log10(expression + 1)",
                  main = paste("Highly Expressed hosts", name), log = "x")
      print(ggDensityTerrain(tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$DistToIntron3,
                  tabs[[name]][tabs[[name]]$MANE_all_exp.norm >= 10, ]$exp.log10, log_x = TRUE,
                  xlab = "Distance to intron 3'", ylab = "Log10(expression + 1)",
                  title = paste("Highly Expressed hosts", name)))
    }
  }

  # Heatscatter: RelIntronPos vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$RelIntronPos, tabs[[name]]$exp.norm, cor = T,
                xlab = "Relative position in the intron",
                ylab = "Normalised expression",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$RelIntronPos, tabs[[name]]$exp.norm,
                xlab = "Relative position in the intron", ylab = "Normalised expression", title = name))
  } else{
    dfSum[name,"CorrExpRelPosIntron"]=cor.test(tabs[[name]]$RelIntronPos, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpRelPosIntron"]=cor.test(tabs[[name]]$RelIntronPos, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpRelPosIntron"]=sum(!is.na(tabs[[name]]$RelIntronPos) & !is.na(tabs[[name]]$exp.log10),na.rm=T)

    heatscatter(tabs[[name]]$RelIntronPos, tabs[[name]]$exp.log10, cor = T,
                xlab = "Relative position in the intron",
                ylab = "Log10(expression + 1)",
                main = name)
    print(ggDensityTerrain(tabs[[name]]$RelIntronPos, tabs[[name]]$exp.log10,
                xlab = "Relative position in the intron", ylab = "Log10(expression + 1)", title = name))
  }

  # Heatscatter: MANE_intron_length vs exp.norm
  if (name %in% names_with_P_list) {
    heatscatter(tabs[[name]]$MANE_intron_length, tabs[[name]]$exp.norm, cor = T,
                xlab = "Intron length",
                ylab = "Normalized expression",
                main = name, log = "x")
    print(ggDensityTerrain(tabs[[name]]$MANE_intron_length, tabs[[name]]$exp.norm, log_x = TRUE,
                xlab = "Intron length", ylab = "Normalized expression", title = name))
  } else {
    dfSum[name,"CorrExpIntronLen"]=cor.test(tabs[[name]]$MANE_intron_length, tabs[[name]]$exp.log10,method="spearman")$estimate
    dfP[name,"CorrExpIntronLen"]=cor.test(tabs[[name]]$MANE_intron_length, tabs[[name]]$exp.log10,method="spearman")$p.value
    dfCounts[name,"CorrExpIntronLen"]=sum(!is.na(tabs[[name]]$MANE_intron_length) & !is.na(tabs[[name]]$exp.log10),na.rm=T)
    heatscatter(tabs[[name]]$MANE_intron_length, tabs[[name]]$exp.log10, cor = T,
                xlab = "Intron length",
                ylab = "Log10(expression + 1)",
                main = name, log = "x")
    print(ggDensityTerrain(tabs[[name]]$MANE_intron_length, tabs[[name]]$exp.log10, log_x = TRUE,
                xlab = "Intron length", ylab = "Log10(expression + 1)", title = name))
  }

  # Heatscatter: MANE_intron_length vs exp.norm for highly expressed hosts
  if (sum(!is.na(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length)) >= 10) {
    if (name %in% names_with_P_list) {
      heatscatter(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length,
                  tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.norm, cor = T,
                  xlab = "Intron length",
                  ylab = "Normalized expression",
                  main = paste("Highly Expressed hosts", name), log = "x")
      print(ggDensityTerrain(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length,
                  tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.norm, log_x = TRUE,
                  xlab = "Intron length", ylab = "Normalized expression",
                  title = paste("Highly Expressed hosts", name)))
    } else {
      dfSum[name,"CorrExpIntronLenHighHosts"]=cor.test(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length,
                                              tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.log10,method="spearman")$estimate
      dfP[name,"CorrExpIntronLenHighHosts"]=cor.test(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length,
                                                       tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.log10,method="spearman")$p.value
      dfCounts[name,"CorrExpIntronLenHighHosts"]=sum(!is.na(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length) &
                                                       !is.na(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.log10),na.rm=T)

      heatscatter(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length,
                  tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.log10, cor = T,
                  xlab = "Intron length",
                  ylab = "Log10(expression + 1)",
                  main = paste("Highly Expressed hosts", name), log = "x")
      print(ggDensityTerrain(tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$MANE_intron_length,
                  tabs[[name]][which(tabs[[name]]$MANE_all_exp.norm >= 10), ]$exp.log10, log_x = TRUE,
                  xlab = "Intron length", ylab = "Log10(expression + 1)",
                  title = paste("Highly Expressed hosts", name)))
    }
  }

  # Boxplot: Expression of integrations by gene
  if (name %in% names_with_P_list) {
    sp = split(tabs[[name]], tabs[[name]]$MANE_all)
    highExp = quantile(tabs[[name]]$exp.norm, 0.9)
    someExp = quantile(tabs[[name]]$exp.norm, 0.75)
    topGenes = sp[sapply(sp, nrow) >= 5]
    if (length(topGenes) >= 2) {
      topGenes = topGenes[order(sapply(topGenes, nrow), decreasing = T)]
      names(topGenes) = trans2gene[match(names(topGenes), trans2gene$V2), "V1"]
      par(mfrow = c(2, 1))
      par(mar = c(5, 5, 5, 5))
      boxplot(lapply(topGenes, "[[", "exp.norm"), las = 2, main = name, ylab = "Integration normalized expression")
      abline(h = c(highExp, someExp), col = "lightgray")
      for (geneInd in 1:length(topGenes)) {
        points(rep(geneInd, nrow(topGenes[[geneInd]])), topGenes[[geneInd]]$exp.norm, pch = 20)
      }

      # Barplot: Host gene expression by gene
      barplot(0.5 + sapply(sapply(topGenes, "[[", "MANE_all_exp"), mean), log = "y", las = 2, main = name, ylab = "Host expression")
      par(mfrow = c(2, 2))
      for (gene in names(topGenes)) {
        plot(topGenes[[gene]]$RelPos, topGenes[[gene]]$exp.norm, main = gene, xlab = "Relative position in the gene", ylab = "Normalized expression", pch = 20)
      }
    } else {
      sp = split(tabs[[name]], tabs[[name]]$MANE_all)
      highExp = quantile(tabs[[name]]$exp.log10, 0.9)
      someExp = quantile(tabs[[name]]$exp.log10, 0.75)
      topGenes = sp[sapply(sp, nrow) >= 5]
      if (length(topGenes) >= 2) {
        topGenes = topGenes[order(sapply(topGenes, nrow), decreasing = T)]
        names(topGenes) = trans2gene[match(names(topGenes), trans2gene$V2), "V1"]
        par(mfrow = c(2, 1))
        par(mar = c(5, 5, 5, 5))
        boxplot(lapply(topGenes, "[[", "exp.log10"), las = 2, main = name, ylab = "Integration log10(expression + 1)")
        abline(h = c(highExp, someExp), col = "lightgray")
        for (geneInd in 1:length(topGenes)) {
          points(rep(geneInd, nrow(topGenes[[geneInd]])), topGenes[[geneInd]]$exp.log10, pch = 20)
        }

        # Barplot: Host gene expression by gene
        barplot(0.5 + sapply(sapply(topGenes, "[[", "MANE_all_exp"), mean), log = "y", las = 2, main = name, ylab = "Host expression")
        par(mfrow = c(2, 2))
        for (gene in names(topGenes)) {
          plot(topGenes[[gene]]$RelPos, topGenes[[gene]]$exp.log10, main = gene, xlab = "Relative position in the gene", ylab = "Log10(expression + 1)", pch = 20)
        }
      }
    }
  }

}
dev.off()
#plot_name="/home/projects/ulitsky/shared/data/human/TRIP/IntronTRIP_run2/plots/test.pdf"
write.table(dfSum,paste0(plot_name,".summary.txt"),sep="\t",quote=F)
write.table(dfCounts,paste0(plot_name,".counts.txt"),sep="\t",quote=F)
write.table(dfP,paste0(plot_name,".pvalue.xt"),sep="\t",quote=F)
pdf(paste0(plot_name,"summary.pdf"), useDingbats = F, width = 14)
summaryPlots = list()
for (name in names(dfSum))
{
  myOrder = order(dfSum[[name]])
  plot_df = data.frame(dataset = factor(niceNames[myOrder], levels = niceNames[myOrder]),
                        value = dfSum[myOrder, name])
  p = ggplot(plot_df, aes(x = dataset, y = value)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = name, x = NULL, y = name) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.5)))
  summaryPlots[[name]] = p
}
# 4 plots per page
for (pageStart in seq(1, length(summaryPlots), by = 4)) {
  pageEnd = min(pageStart + 3, length(summaryPlots))
  print(ggarrange(plotlist = summaryPlots[pageStart:pageEnd], ncol = 2, nrow = 2))
}
pheatmap(na.omit(dfSum))
pheatmap(dfSum,cluster_rows = F,cluster_cols=F)
dev.off()
