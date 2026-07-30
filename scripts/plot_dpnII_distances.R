# =============================================================================
# plot_dpnII_distances.R
# =============================================================================
# Compares the real per-barcode "distance to nearest downstream DpnII site"
# (written by stage 06, scripts/06_get_genuine_barcodes.R, as
# <mapName>.txt.retained_barcodes_genuine.dnpII_dists.txt) against a random-genomic
# background: nSample TTAA sites sampled at random from the real hg38 TTAA reference,
# each assigned a random +/- strand (piggyBac has no inherent orientation preference),
# with the identical nearest-downstream-DpnII-site distance computed for each -- the
# same "what would this look like at a random insertable site" background used
# elsewhere in this pipeline for TTAA enrichment checks.
#
# Usage
#   Rscript plot_dpnII_distances.R <path> <distFiles> <ttaaBedPath> <dpniiBedPath> <gapsBedPath> <outPdf> [nSample]
#
#   path         directory containing both the *.dnpII_dists.txt file(s) and where outPdf
#                is written (unless outPdf is an absolute path)
#   distFiles    comma-separated list of *.dnpII_dists.txt filenames (relative to path) --
#                e.g. a single mapName's file, or several pooled into one histogram
#   ttaaBedPath  hg38 TTAA reference bed (same format/file as stage 04's TTAA filter)
#   dpniiBedPath hg38 DpnII digest-fragment bed (same file stage 06 reads)
#   gapsBedPath  hg38 assembly-gaps bed (chrom, start, end -- N stretches/unassembled
#                regions), used to drop any (site, downstream DpnII site) pair with a gap
#                physically between them, since the "nearest downstream DpnII site" the
#                reference reports for those is not meaningful (the true nearest site, if
#                any, could be on the far side of unknown/unassembled sequence)
#   outPdf       output PDF path
#   nSample      number of random TTAA sites to sample for the background (default 10000)
# =============================================================================

path = "/home/projects/ulitsky/shared/data/human/TRIP/single_combined/"
distFiles = "5ss_noP.txt.retained_barcodes_genuine.dnpII_dists.txt"
ttaaBedPath = "/home/projects/ulitsky/igoru/data/human/hg38_ttaa.bed"
dpniiBedPath = "/home/projects/ulitsky/igoru/data/human/hg38_dpnii.bed"
gapsBedPath = "/home/projects/ulitsky/igoru/data/human/hg38.gaps.bed.gz"
outPdf = "plots/dpnII_distances.pdf"
nSample = 10000

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs) > 0) {
  path = myArgs[1]
  distFiles = myArgs[2]
  ttaaBedPath = myArgs[3]
  dpniiBedPath = myArgs[4]
  gapsBedPath = myArgs[5]
  outPdf = myArgs[6]
  if (length(myArgs) >= 7) nSample = as.integer(myArgs[7])
}

set.seed(1)

# Same DpnII-site extraction and nearest-downstream-distance logic as stage 06
# (self-contained here rather than sourced, matching this codebase's convention of
# not sharing code across per-script files).
print(paste("Loading DpnII sites from", dpniiBedPath))
dpnii = read.table(dpniiBedPath, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = c("character","integer","NULL","NULL","NULL","NULL"))
colnames(dpnii) = c("chrom","start")
dpnii = dpnii[order(dpnii$chrom, dpnii$start), ]
dpnii = dpnii[duplicated(dpnii$chrom), ]                # drop each chrom's first (boundary) fragment
dpniiByChrom = lapply(split(dpnii$start, dpnii$chrom), sort)
print(paste("Loaded", nrow(dpnii), "DpnII sites on", length(dpniiByChrom), "contigs"))

nearestDownstreamDpnII <- function(chrom, pos, strand, byChrom)
{
  vec = byChrom[[chrom]]
  if (is.null(vec)) return(NA_real_)
  if (strand == "+") {
    idx = findInterval(pos, vec) + 1L
    if (idx > length(vec)) return(NA_real_)
    vec[idx] - pos
  } else {
    idx = findInterval(pos - 1L, vec)
    if (idx < 1L) return(NA_real_)
    pos - vec[idx]
  }
}

# Assembly gaps (N stretches / unassembled sequence): a (site, downstream DpnII site) pair
# whose interval spans a real gap is dropped rather than trusted, regardless of its
# distance -- the reference genuinely has no sequence there, so "nearest downstream DpnII
# site" is not a meaningful measurement for that pair.
print(paste("Loading assembly gaps from", gapsBedPath))
gapsCon = if (grepl("\\.gz$", gapsBedPath)) gzfile(gapsBedPath) else gapsBedPath
gaps = read.table(gapsCon, sep = "\t", stringsAsFactors = FALSE)[, 1:3]
colnames(gaps) = c("chrom","start","end")
gapsByChrom = split(gaps[, c("start","end")], gaps$chrom)
print(paste("Loaded", nrow(gaps), "assembly gaps on", length(gapsByChrom), "contigs"))

hasGapBetween <- function(chrom, pos, strand, dist, byChrom)
{
  if (is.na(dist)) return(FALSE)
  g = byChrom[[chrom]]
  if (is.null(g)) return(FALSE)
  if (strand == "+") { segStart = pos; segEnd = pos + dist }
  else { segStart = pos - dist; segEnd = pos }
  any(g$start < segEnd & g$end > segStart)
}

# A handful of distances come out implausibly large (multi-kb to tens of kb) -- almost
# certainly assembly gaps (N stretches) or other low-complexity/repeat-masked regions
# where the reference genuinely has no GATC for a long stretch, not a real biological
# signal. Cap at 10x the by-chance expected distance (assuming a gap-free genome at
# genomeGC% GC content, independent bases: P(GATC) = (GC/2)^2 * ((1-GC)/2)^2, expected
# distance = 1/P(GATC)) and drop anything beyond that rather than let a few genome-gap
# outliers distort the histograms/means. Applied to both real and background distances
# identically -- the real dnpII_dists.txt files already have this applied by stage 06,
# but re-applying here is a no-op on those and correctly still applies to the freshly
# computed random-TTAA background below.
genomeGC = 0.44
pGATC = (genomeGC/2)^2 * ((1-genomeGC)/2)^2
expectedDist = 1/pGATC
maxPlausibleDist = 10 * expectedDist
print(paste0("Expected by-chance DpnII spacing at ", genomeGC*100, "% GC (no Ns): ",
             round(expectedDist,1), " bp; dropping distances beyond 10x that (",
             round(maxPlausibleDist,1), " bp)"))

# --- Real distances: pool every listed dnpII_dists.txt file ---
print("Loading real barcode DpnII distances")
realTabs = lapply(strsplit(distFiles, ",")[[1]], function(f) read.table(paste0(path, f), sep = "\t", header = TRUE, stringsAsFactors = FALSE))
realAll = do.call(rbind, realTabs)
realGapped = mapply(hasGapBetween, realAll$chr, realAll$start, realAll$strand, realAll$dpnII_dist,
                     MoreArgs = list(byChrom = gapsByChrom))
print(paste(sum(realGapped, na.rm = TRUE), "real distances dropped for spanning an assembly gap"))
realDist = realAll$dpnII_dist
realDist = realDist[!is.na(realDist) & realDist <= maxPlausibleDist & !realGapped]
print(paste(length(realDist), "real genuine-barcode distances loaded"))

# --- Background: nSample random TTAA sites, random strand, same distance metric ---
print(paste("Loading TTAA reference from", ttaaBedPath))
ttaa = read.table(ttaaBedPath, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = c("character","integer","NULL","NULL","NULL","NULL"))
colnames(ttaa) = c("chrom","start")
print(paste("Sampling", nSample, "random TTAA sites (of", nrow(ttaa), ") with a random strand each"))
sampleIdx = sample(nrow(ttaa), min(nSample, nrow(ttaa)))
randChrom = ttaa$chrom[sampleIdx]
randPos = ttaa$start[sampleIdx]
randStrand = sample(c("+","-"), length(sampleIdx), replace = TRUE)
randDist = mapply(nearestDownstreamDpnII, randChrom, randPos, randStrand, MoreArgs = list(byChrom = dpniiByChrom))
randGapped = mapply(hasGapBetween, randChrom, randPos, randStrand, randDist, MoreArgs = list(byChrom = gapsByChrom))
print(paste(sum(randGapped, na.rm = TRUE), "random-TTAA distances dropped for spanning an assembly gap"))
randDist = randDist[!is.na(randDist) & randDist <= maxPlausibleDist & !randGapped]
print(paste(length(randDist), "random-TTAA background distances computed"))

print(paste("Real: median", median(realDist), " mean", round(mean(realDist),1)))
print(paste("Random TTAA background: median", median(randDist), " mean", round(mean(randDist),1)))

dir.create(dirname(paste0(path, outPdf)), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(path, outPdf), width = 10, height = 5.5)
par(mfrow = c(1,2), oma = c(0,0,2.5,0))

capLabel = paste0("outlier-exclusion cap = ", round(maxPlausibleDist,1),
                   " bp (10x the ", round(expectedDist,1), " bp expected-by-chance spacing at ",
                   genomeGC*100, "% GC)")

# Linear scale, zoomed to the 99th percentile of the pooled data (a SEPARATE, purely
# cosmetic zoom -- not the outlier-exclusion cap above) so a handful of distant points
# don't compress the whole histogram into one bar. This zoom is often well below
# maxPlausibleDist, so values beyond the visible x-axis here are real, retained data, not
# a sign the outlier cap failed -- see the log10 panel for the full range, and capLabel
# above for the actual exclusion threshold.
capAt = quantile(c(realDist, randDist), 0.99)
breaksLinear = seq(0, capAt, length.out = 41)
histReal = hist(pmin(realDist, capAt), breaks = breaksLinear, plot = FALSE)
histRand = hist(pmin(randDist, capAt), breaks = breaksLinear, plot = FALSE)
yMax = max(histReal$density, histRand$density)
plot(histRand, freq = FALSE, col = rgb(0.5,0.5,0.5,0.5), border = NA, ylim = c(0, yMax),
     xlab = "Distance to nearest downstream DpnII site (bp)",
     main = paste0("Linear scale (zoomed to pooled p99 = ", round(capAt), " bp)"))
plot(histReal, freq = FALSE, col = rgb(0.85,0.2,0.2,0.5), border = NA, add = TRUE)
if (maxPlausibleDist <= capAt) {
  abline(v = maxPlausibleDist, lty = 2)
  text(maxPlausibleDist, yMax, "cap", pos = 2, cex = 0.8)
}
legend("topright", legend = c(paste0("Real (n=", length(realDist), ")"),
                               paste0("Random TTAA (n=", length(randDist), ")")),
       fill = c(rgb(0.85,0.2,0.2,0.5), rgb(0.5,0.5,0.5,0.5)), border = NA, bty = "n")

# log10(dist+1), full range (up to maxPlausibleDist -- everything beyond that was already
# excluded above), so the outlier cap itself is visible as a vertical line here.
logReal = log10(realDist + 1)
logRand = log10(randDist + 1)
breaksLog = seq(0, max(c(logReal, logRand)), length.out = 41)
histRealLog = hist(logReal, breaks = breaksLog, plot = FALSE)
histRandLog = hist(logRand, breaks = breaksLog, plot = FALSE)
yMaxLog = max(histRealLog$density, histRandLog$density)
plot(histRandLog, freq = FALSE, col = rgb(0.5,0.5,0.5,0.5), border = NA, ylim = c(0, yMaxLog),
     xlab = "log10(distance + 1) (bp)", main = "Log scale (full retained range)")
plot(histRealLog, freq = FALSE, col = rgb(0.85,0.2,0.2,0.5), border = NA, add = TRUE)
abline(v = log10(maxPlausibleDist + 1), lty = 2)
text(log10(maxPlausibleDist + 1), yMaxLog, "cap", pos = 2, cex = 0.8)
legend("topright", legend = c(paste0("Real (n=", length(realDist), ")"),
                               paste0("Random TTAA (n=", length(randDist), ")")),
       fill = c(rgb(0.85,0.2,0.2,0.5), rgb(0.5,0.5,0.5,0.5)), border = NA, bty = "n")

mtext(capLabel, outer = TRUE, cex = 0.9, line = 0.5)

dev.off()
print(paste("Wrote", paste0(path, outPdf)))
