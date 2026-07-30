# =============================================================================
# Stage 14 — combine spliced/unspliced/total combined_counts.txt triplets into
# splicing-ratio readouts
# =============================================================================
# Faithful port of igoru's ~/R/trip/IntronTrip/combineSplicedUnsplicedPairs_uniform.R --
# no logic changes. Runs after stage 07: each row of ### SPLICED_UNSPLICED (written out as
# spliced_unspliced.txt by run_pipeline.sh) names three ALREADY-PRODUCED stage 07 outputs
# (spliced/unspliced/total <expName>.combined_counts.txt, same barcode-position rows since
# they share the same map/norm) and an output name prefix; this stage derives three ratio
# tables per row.
#
# Usage
#   Rscript 14_combineSplicedUnsplicedPairs_uniform.R <path>
#
# Input   <path>/spliced_unspliced.txt -- outName<TAB>splicedExpName<TAB>unsplicedExpName<TAB>totalExpName
#         <path>/<splicedExpName|unsplicedExpName|totalExpName>.combined_counts.txt (stage 07)
# Outputs per row, all under <path>/<outName>.*:
#   .spliced_unspliced.combined_counts.txt            exp.log10 = spliced/(spliced+unspliced), NA if <20 reads
#   .spliced_to_total.combined_counts.txt              rows with spliced+total>=10 kept; exp.log10 = spliced/(spliced+total), NA if <20
#   .spliced_and_unspliced_to_total.combined_counts.txt  further filtered to sum>=10; exp.log10 = (spliced+unspliced)/(spliced+unspliced+total)
# =============================================================================

path = "/home/projects/ulitsky/shared/data/human/TRIP/SON_combined/"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path=myArgs[1]
}

pairs = read.table(paste0(path,"spliced_unspliced.txt"),sep="\t")

for (ind in 1:nrow(pairs))
{
  print(ind)
  spliced=read.table(paste0(path,pairs[ind,2],".combined_counts.txt"),sep="\t",header=T)
  print(nrow(spliced))
  unspliced=read.table(paste0(path,pairs[ind,3],".combined_counts.txt"),sep="\t",header=T)
  print(nrow(unspliced))
  total=read.table(paste0(path,pairs[ind,4],".combined_counts.txt"),sep="\t",header=T)
  print(nrow(total))

  stopifnot(nrow(spliced)==nrow(unspliced),nrow(spliced)==nrow(total))
  merged = spliced

  merged$exp_unspliced = unspliced$exp
  merged$exp_total = total$exp
  merged$exp.log10=(merged$exp)/(merged$exp+merged$exp_unspliced)
  merged$exp.log10[merged$exp+merged$exp_unspliced<20]=NA
  write.table(merged,paste0(path,pairs[ind,1],".spliced_unspliced.combined_counts.txt"),row.names = T,col.names = T,sep="\t",quote=F)

  merged=merged[merged$exp+merged$exp_total>=10,]
  merged$exp.log10=(merged$exp)/(merged$exp+merged$exp_total)
  merged$exp.log10[merged$exp+merged$exp_total<20]=NA
  write.table(merged,paste0(path,pairs[ind,1],".spliced_to_total.combined_counts.txt"),row.names = T,col.names = T,sep="\t",quote=F)

  merged=merged[merged$exp+merged$exp_unspliced+merged$exp_total>=10,]
  merged$exp.log10=(merged$exp_unspliced+merged$exp)/(merged$exp+merged$exp_unspliced+merged$exp_total)
  write.table(merged,paste0(path,pairs[ind,1],".spliced_and_unspliced_to_total.combined_counts.txt"),row.names = T,col.names = T,sep="\t",quote=F)
}
