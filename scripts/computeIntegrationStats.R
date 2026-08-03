## Usage: Rscript computeIntegrationStats.R XXX
## Reads XXX.annotated_table_used.txt (stage 10's per-construct annotated table, restricted
## to barcodes with >=2 total normalization reads across every COMBOS normSample for that
## construct -- one row per TRIP integration/barcode; ragged header -> first data column is
## the barcode/rownames, as used throughout this pipeline -- see readAnnotTable() below) and
## reports:
##  - average integration frequency (mean/median distance between consecutive integrations
##    on the same chromosome)
##  - how many MANE transcripts and how many GENCODE genes have >=1 integration (out of the
##    total number that exist)
##  - what fraction of integrations fall in a) MANE introns, b) MANE exons, c) GENCODE exons
##    not already counted in a/b, d) GENCODE introns not already counted in a-c, e) antisense to
##    a MANE transcript and not already counted in a-d (mutually exclusive, in that priority
##    order; everything left over is "other/intergenic")
##
## a) and b) use the MANE_intron/MANE_exonic columns already present in the annotated table
## file (sense-strand transcript/intron overlap, computed by the annotation pipeline itself).
## c) and d) are computed fresh here via a sense-strand interval join against the full GENCODE
## v42 exon/intron annotation (gencode.v42.annotation.exons.bed / .introns.bed), which covers
## every GENCODE transcript including the MANE-select ones -- so the union of exon+intron hits
## there is also what's used as the "total GENCODE genes hit" denominator, giving a single
## self-consistent gene set rather than mixing MANE's RefSeq IDs with GENCODE's Ensembl IDs.
##
## This is a standalone visualization/reporting script, not a numbered pipeline stage -- run
## it by hand against a construct's stage 10 output whenever you want the report/pie chart,
## same convention as plot_dpnII_distances.R and plotSplicingEfficiencies_uniform.R.

suppressMessages(library(data.table))

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript computeIntegrationStats.R XXX  (reads XXX.annotated_table_used.txt)")
xxx = args[1]
infile = paste0(xxx, ".annotated_table_used.txt")

maneBedFile = "/home/projects/ulitsky/shared/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed"
gencodeExonsFile = "/home/projects/ulitsky/shared/data/human/gencode.v42.annotation.exons.bed"
gencodeIntronsFile = "/home/projects/ulitsky/shared/data/human/gencode.v42.annotation.introns.bed"
geneInfoFile = "/home/projects/ulitsky/shared/data/human/gencode.v42.gene_info.txt"

# Ragged-header convention used throughout this pipeline: the header line has one fewer
# field than the data rows (the extra leading data field is the barcode -> rownames).
readAnnotTable <- function(file) {
  headerLine = strsplit(readLines(file, n = 1), "\t")[[1]]
  x = fread(file, sep = "\t", header = FALSE, skip = 1, data.table = FALSE)
  rownames(x) = x[[1]]
  x = x[, -1, drop = FALSE]
  names(x) = headerLine
  x
}

cat("Reading", infile, "...\n")
tab = readAnnotTable(infile)
n = nrow(tab)
nUnique = length(unique(paste(tab$Chr, tab$Pos, tab$Strand)))
cat(format(n, big.mark = ","), "total integrations (barcodes)\n")
cat(format(nUnique, big.mark = ","), "unique integrations (distinct Chr:Pos:Strand)\n\n")

## ---- a) integration frequency: mean/median distance between consecutive integrations,
## computed within each chromosome (never across a chromosome boundary) ----
gaps = unlist(lapply(split(tab$Pos, tab$Chr), function(p) diff(sort(p))))
cat("=== Integration frequency ===\n")
cat("Mean distance between consecutive integrations (same chromosome):", format(mean(gaps), big.mark = ","), "bp\n")
cat("Median distance between consecutive integrations (same chromosome):", format(median(gaps), big.mark = ","), "bp\n\n")

## ---- b) MANE transcripts / GENCODE genes with >=1 integration ----
maneBed = fread(maneBedFile, sep = "\t", header = FALSE, data.table = FALSE)
maneTotalTranscripts = length(unique(maneBed$V4))
maneHitTranscripts = length(unique(na.omit(tab$MANE_all)))

## ---- classify each integration: MANE intron > MANE exon > GENCODE exon > GENCODE intron >
## antisense-to-MANE > other ----
isManeIntron = !is.na(tab$MANE_intron)
isManeExon = !is.na(tab$MANE_exonic) & !isManeIntron

# Union of all three antisense columns -- MANE_antisense (whole transcript body) is NOT a
# perfect superset of MANE_intron_antisense/MANE_exonic_antisense (confirmed on real data:
# ~0.4% of intron-antisense-flagged rows have MANE_antisense == NA), so any one of the three
# being non-NA counts as "antisense to a MANE transcript".
isAntisenseToMane = !is.na(tab$MANE_antisense) | !is.na(tab$MANE_exonic_antisense) | !is.na(tab$MANE_intron_antisense)

# Sense-strand interval join against the full GENCODE exon/intron annotation, restricted to
# rows not already classified as MANE intron/exon.
toCheck = which(!isManeIntron & !isManeExon)
query = data.table(Chr = tab$Chr[toCheck], Strand = tab$Strand[toCheck],
                    start = tab$Pos[toCheck] - 1L, end = tab$Pos[toCheck],
                    qidx = toCheck)
setkey(query, Chr, Strand, start, end)

loadGencodeBed <- function(file) {
  x = fread(file, sep = "\t", header = FALSE, select = 1:6,
            col.names = c("Chr", "start", "end", "id", "score", "Strand"), data.table = TRUE)
  x$geneTranscript = sub("\\..*", "", x$id)
  setkey(x, Chr, Strand, start, end)
  x
}

cat("Loading GENCODE exon/intron annotation...\n")
gencodeExons = loadGencodeBed(gencodeExonsFile)
gencodeIntrons = loadGencodeBed(gencodeIntronsFile)

exonHits = foverlaps(query, gencodeExons, nomatch = NULL)
intronHits = foverlaps(query, gencodeIntrons, nomatch = NULL)

isGencodeExon = rep(FALSE, n)
isGencodeExon[unique(exonHits$qidx)] = TRUE

isGencodeIntron = rep(FALSE, n)
isGencodeIntron[unique(intronHits$qidx)] = TRUE
isGencodeIntron = isGencodeIntron & !isGencodeExon & !isManeIntron & !isManeExon
isGencodeExon = isGencodeExon & !isManeIntron & !isManeExon

# Only assign the antisense-to-MANE category from what's left after the first four.
isAntisenseToMane = isAntisenseToMane & !isManeIntron & !isManeExon & !isGencodeExon & !isGencodeIntron

## GENCODE genes hit by ANY integration (any category), via ENST -> ENSG on the same
## exon/intron overlap (not restricted to c/d only, since MANE-select transcripts are also
## present in this annotation under their ENST id -- see script header comment).
geneInfo = fread(geneInfoFile, sep = "\t", header = TRUE, data.table = FALSE)
geneInfo = geneInfo[!duplicated(geneInfo[["Transcript stable ID"]]), ]
rownames(geneInfo) = geneInfo[["Transcript stable ID"]]

allHitTranscripts = unique(c(exonHits$geneTranscript, intronHits$geneTranscript))
allHitGenes = unique(na.omit(geneInfo[allHitTranscripts, "Gene stable ID"]))
gencodeTotalGenes = length(unique(geneInfo[["Gene stable ID"]]))

cat("\n=== Gene/transcript coverage ===\n")
cat("MANE transcripts with >=1 integration:", maneHitTranscripts, "/", maneTotalTranscripts,
    sprintf("(%.2f%%)\n", 100 * maneHitTranscripts / maneTotalTranscripts))
cat("GENCODE genes with >=1 integration:", length(allHitGenes), "/", gencodeTotalGenes,
    sprintf("(%.2f%%)\n\n", 100 * length(allHitGenes) / gencodeTotalGenes))

## ---- c) mutually-exclusive breakdown ----
categories = ifelse(isManeIntron, "MANE intron",
              ifelse(isManeExon, "MANE exon",
              ifelse(isGencodeExon, "GENCODE exon (non-MANE)",
              ifelse(isGencodeIntron, "GENCODE intron (non-MANE)",
              ifelse(isAntisenseToMane, "Antisense to MANE transcript", "other/intergenic")))))
tabCounts = table(factor(categories, levels = c("MANE intron", "MANE exon",
                                                 "GENCODE exon (non-MANE)",
                                                 "GENCODE intron (non-MANE)",
                                                 "Antisense to MANE transcript",
                                                 "other/intergenic")))

cat("=== Integration location breakdown ===\n")
for (cat_name in names(tabCounts)) {
  cnt = tabCounts[[cat_name]]
  cat(sprintf("%-28s %8s  (%.2f%%)\n", cat_name, format(cnt, big.mark = ","), 100 * cnt / n))
}

## ---- pie chart of the breakdown ----
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
suppressMessages(library(ggrepel))

catLevels = c("MANE intron", "MANE exon", "GENCODE exon (non-MANE)",
              "GENCODE intron (non-MANE)", "Antisense to MANE transcript", "other/intergenic")
# Fixed categorical hue order (blue, orange, aqua, yellow, magenta, green) -- assigned by
# category identity, not by slice size, and never reused/recycled across runs.
catColors = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")

pieDF = data.frame(category = factor(catLevels, levels = catLevels),
                    count = as.numeric(tabCounts[catLevels]))
pieDF$pct = 100 * pieDF$count / sum(pieDF$count)
# Every slice gets a label (no size threshold) -- category name, %, and n, in place of a
# separate legend.
pieDF$label = sprintf("%s\n%.1f%%  (n=%s)", pieDF$category, pieDF$pct, format(pieDF$count, big.mark = ","))

# Built with ggforce::geom_arc_bar (ordinary Cartesian wedge polygons, angle 0 = 12
# o'clock, increasing clockwise) rather than coord_polar()+geom_col, so ggrepel's usual
# Cartesian-space repelling/leader-line logic applies to the labels -- coord_polar's
# transformed space doesn't work with ggrepel. Every label sits just outside the pie and
# is nudged apart from its neighbors automatically, which matters most for small
# neighboring slices (e.g. "MANE exon"/"GENCODE exon (non-MANE)" here) that would
# otherwise overlap if drawn inside their own sliver.
pieDF$fraction = pieDF$count / sum(pieDF$count)
pieDF$end = cumsum(pieDF$fraction) * 2 * pi
pieDF$start = c(0, head(pieDF$end, -1))
pieDF$mid = (pieDF$start + pieDF$end) / 2
r = 1
# ggrepel draws each leader line from the exact (x,y) data point to the label, so the
# anchor has to sit ON the pie's own edge (r itself) for the line to actually touch the
# pie -- an anchor placed further out (as a previous version of this script did, at 1.5x
# r) leaves a floating gap between the pie and where the line starts. Labels are instead
# pushed outward via repulsion strength (force/force_pull/box.padding below), not by
# moving the anchor.
pieDF$labelX = r * sin(pieDF$mid)
pieDF$labelY = r * cos(pieDF$mid)

p = ggplot(pieDF) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = r, start = start, end = end, fill = category),
               color = "#fcfcfb", linewidth = 0.8) +
  geom_label_repel(aes(x = labelX, y = labelY, label = label),
                    size = 3, lineheight = 0.9, fontface = "bold",
                    segment.color = "grey40", min.segment.length = 0,
                    label.size = NA, fill = scales::alpha("white", 0.85),
                    box.padding = 0.6, point.padding = 0.3,
                    force = 4, force_pull = 0.2) +
  scale_fill_manual(values = catColors, breaks = catLevels, guide = "none") +
  coord_fixed(clip = "off") +
  expand_limits(x = c(-2.4, 2.4), y = c(-2.4, 2.4)) +
  labs(title = "TRIP integration locations",
       subtitle = paste0(basename(xxx), "  (n = ", format(n, big.mark = ","), ")")) +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "#52514e"))

dir.create("plots", showWarnings = FALSE)
pieFile = paste0("plots/", xxx, ".integration_categories_pie.pdf")
ggsave(pieFile, p, width = 7, height = 5.5)
cat("\nPie chart saved to", pieFile, "\n")
