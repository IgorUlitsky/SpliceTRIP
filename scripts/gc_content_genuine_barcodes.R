# =============================================================================
# gc_content_genuine_barcodes — per-construct GC-content vs expression report
# =============================================================================
# Visualization/QC script, not a numbered pipeline stage -- run by hand against stage
# 07/10's output whenever you want this report, same convention as
# plotSplicingEfficiencies_uniform.R and computeIntegrationStats.R.
#
# Usage
#   Rscript gc_content_genuine_barcodes.R <path> <mapName>.txt <namesExp_file> <outFile.pdf>
#
# Inputs   <mapName>.txt.annotated_tables.txt   (stage 10)
#          <expName>.combined_counts.txt        (stage 07, or any pipeline producing the
#                                                same columns -- e.g. a splicing-ratio
#                                                exp.log10, auto-detected and preserved;
#                                                see examples/exp_lists/ for real ones)
#   namesExp_file: one <expName> per line, same format as annotateTRIP_positions_genuine_
#                  barcodes.R/annotateTRIPGeneBodies_genuine_barcodes.R's.
#
# Outputs  <outFile.pdf>   per-library: GC% density histogram, expression-vs-GC% scatter
#                          (colored by Low/Average/High GC category), and expression-by-
#                          GC-category boxplot (pairwise t-tests between the 3 categories)
#
# Run once PER construct (loop over every distinct ### FASTQ_MAPS mapName), same as
# annotateTRIP_positions_genuine_barcodes.R/annotateTRIPGeneBodies_genuine_barcodes.R --
# unlike comparePairwiseSharedPositions_uniform.R, which takes the whole construct list in
# one call.
#
# This is a verbatim port of igoru's `gc_content_genuine_barcodes.R`
# (/home/projects/ulitsky/lidiag/R/scripts/decent_scripts/), aside from an added doc
# header and one real bug fix (applied to the original script directly too): the
# unconditional `exp.log10 = log10(1+exp)` was silently overwriting the pre-existing,
# differently-scaled exp.log10 that *.spliced_unspliced.combined_counts.txt files already
# carry (a [0,1]-bounded splicing ratio) -- confirmed run in production against exactly
# such files (real gc_content_IntronTRIP_spliced_unspliced.pdf output on disk). Same fix
# already applied to stage 07/10's and annotateTRIP_positions_genuine_barcodes.R's/
# annotateTRIPGeneBodies_genuine_barcodes.R's equivalent merges. No hardcoded external
# paths to extract -- like those two scripts, it only reads other pipeline stages' own
# output.
# =============================================================================
library(ggplot2)
library(ggpubr)
library("RColorBrewer")

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined/"
namesExp<-read.table(paste0(path,"names.5ss.txt"))$V1


rep_with_promoter <- c('5ss_Expression_rep1_S56 5ss_Expression_rep2_S57',
                       'CPA_Expression_rep2_S54 CPA_Expression_rep3_S55')


merged = list()

# Arguments from the terminal
myArgs = commandArgs(trailingOnly = T)
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

names_with_P_list = c('5ss_Expression_rep1_S56', '5ss_Expression_rep2_S57',
                      'CPA_Expression_rep2_S54', 'CPA_Expression_rep3_S55')


pdf(plot_name, useDingbats = FALSE, width = 12)

annotTable =  read.table(paste0(path, nameMap, ".annotated_tables.txt"), sep = "\t", quote = "")

for (name in namesExp) {
  print(name)
  countsTab =  read.table(paste0(path, name, ".combined_counts.txt"), sep = "\t", quote = "")
  countsTab$Barcode = rownames(countsTab)
  annotTable$Barcode = rownames(annotTable)
  myMerged = merge(annotTable,countsTab,by=c("Barcode","Chr","Pos","Strand"))
  # Some combined_counts.txt files (e.g. *.spliced_unspliced.combined_counts.txt) already
  # carry their own exp.log10 -- a splicing ratio bounded [0,1], NA when undefined -- which
  # is a completely different metric from log10(1+exp) below; don't silently overwrite it.
  if (!("exp.log10" %in% names(myMerged))) {
    myMerged$exp.log10 = log(1+myMerged$exp,10)
  }
  myMerged$exp.norm = log((1 + myMerged$exp) / (1 + myMerged$norm), 2)
  rownames(myMerged)=myMerged$Barcode
  myMerged$Row.names=NULL


  merged[[name]] <- myMerged
  #Density
  if (name %in% names_with_P_list) {
    p1 <- ggplot(merged[[name]], aes(x = GC)) +
      geom_histogram(aes(y = ..density..), binwidth = 1, fill = "blue", color = "black", alpha = 0.5) +
      geom_density(color = "red", size = 1) +
      labs(title = paste("Distribution of GC composition values +/- 500 bp from the insertion site,", name, sep=' '),
           x = "GC%",
           y = "Density") +
      theme_minimal()
  } else {
    p1 <-ggplot(merged[[name]], aes(x = GC)) +
      geom_histogram(aes(y = ..density..), binwidth = 1, fill = "blue", color = "black", alpha = 0.5) +
      geom_density(color = "red", size = 1) +
      labs(title = paste("Distribution of GC composition values +/- 500 bp from the insertion site,", name, sep=' '),
           x = "GC%",
           y = "Density") +
      theme_minimal()
  }
  print(p1)

  # Exp ~ GC
  gc_counts <- table(merged[[name]]$GC_category)
  gc_labels <- paste(names(gc_counts), "(", gc_counts, ")", sep = "")
  if (name %in% names_with_P_list) {
    p2 <- ggplot(merged[[name]], aes(x = exp.norm, y = GC, color = GC_category)) +
      geom_point(alpha = 0.6) +
      labs(title = paste("Expression ~ GC%", name, sep=' '),
           x = "exp.norm",
           y = "GC%") +
      theme_minimal() +
      scale_color_manual(values = c("Low GC" = "#0b3954", "Average GC" = "#bfd7ea", "High GC" = "#ff6663"),
                         labels = gc_labels)
  } else {
    p2 <- ggplot(merged[[name]], aes(x = exp.log10, y = GC, color = GC_category)) +
      geom_point(alpha = 0.6) +
      labs(title = paste("Expression ~ GC%", name, sep=' '),
           x = "exp.log10",
           y = "GC%") +
      theme_minimal() +
      scale_color_manual(values = c("Low GC" = "#0b3954", "Average GC" = "#bfd7ea", "High GC" = "#ff6663"),
                         labels = gc_labels)
  }
  print(p2)

  # Exp ~ GC, Boxplot
  if (name %in% names_with_P_list) {
    p3 <- ggplot(merged[[name]], aes(x = GC_category, y = exp.norm, fill = GC_category)) +
      geom_boxplot(alpha = 0.6) +
      labs(title = paste("Boxplot of expression values by GC category", name, sep=' '),
           x = "GC category",
           y = "exp.norm")  +
      theme_minimal()+
      scale_fill_manual(values = c("Low GC" = "#0b3954", "Average GC" = "#bfd7ea", "High GC" = "#ff6663"),
                        labels = gc_labels) +
      stat_compare_means(method = "t.test", label = "p.format",
                         comparisons = list(c("Low GC", "Average GC"),
                                            c("Low GC", "High GC"),
                                            c("Average GC", "High GC")),
                         size = 4)
  } else {
    p3 <- ggplot(merged[[name]], aes(x = GC_category, y = exp.log10, fill = GC_category)) +
      geom_boxplot(alpha = 0.6) +
      labs(title = paste("Boxplot of expression values by GC category", name, sep=' '),
           x = "GC category",
           y = "exp.log10") +
      theme_minimal()+
      scale_fill_manual(values = c("Low GC" = "#0b3954", "Average GC" = "#bfd7ea", "High GC" = "#ff6663"),
                        label = gc_labels) +
      stat_compare_means(method = "t.test", label = "p.format",
                         comparisons = list(c("Low GC", "Average GC"),
                                            c("Low GC", "High GC"),
                                            c("Average GC", "High GC")),
                         size = 4)
  }
  print(p3)
}

dev.off()
