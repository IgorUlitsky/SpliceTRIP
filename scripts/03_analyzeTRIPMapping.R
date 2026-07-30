library(DescTools)
library(ggplot2)
library(cowplot)
library(gridExtra)

#library(Rncc)

# Take in parameters
path = "/home/projects/ulitsky/lidiag/TRIP_data/run2/"
path_to_mapping = "/home/projects/ulitsky/lidiag/TRIP_data/run2/02_mapping/"
name= "5ss_Map_rep1_S64"
minAlignLen = 20
maxAlignLen = 300
barcodeLen = 18
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
    path=myArgs[1]
    path_to_mapping=myArgs[2]
    name=myArgs[3]
    # 4th arg is optional; barcodeLen keeps its default of 18 above when omitted (barcodes
    # in this pipeline are 18nt for most constructs, but not all -- some use 16nt).
    if (length(myArgs)>=4) barcodeLen = as.integer(myArgs[4])
}
barcodes = list()
print(name)
print("Loading barcodes")
tabBarcodes<-read.table(gzfile(paste0(path,name,".barcodes.txt.gz")),sep="\t")
print("Loading positions")
tabStarts<-read.table(gzfile(paste0(path_to_mapping,name,".second_read.bed.gz")),sep="\t")
tabStarts[1,]


print(paste(nrow(tabStarts),"input alignments"))
tabStarts$AlignLen =tabStarts$V3-tabStarts$V2

# ggplot(tabStarts, aes(x = factor(1), y = AlignLen)) +
#   geom_boxplot(fill = "lightblue", color = "black", width = 0.2, position = position_nudge(x = 0.45)) +
#   geom_violin(fill = "lightgreen", color = NA, width = 0.5) +
#   projects(title = name, y = "AlignLen", x = "") +
#   coord_flip() +
#   theme_minimal()

print(paste(sum(tabStarts$AlignLen<minAlignLen)," have alignments shorter than",minAlignLen))
print(paste(sum(tabStarts$AlignLen>maxAlignLen)," have alignments longer than",maxAlignLen))
tabStarts=subset(tabStarts,tabStarts$AlignLen>=minAlignLen & tabStarts$AlignLen<=maxAlignLen)
print(paste(nrow(tabStarts),"valid length alignments"))

# Fixed for 2nd read
# Minus-strand position uses V3 (not V3-1). A prior version of this line used V3-1 on the
# theory that BED end is exclusive so V3-1 is the "true last aligned base" -- reverted after
# it turned out to be wrong. Verified directly against real, already-computed
# .barcodes.bed files (run3/5ss_Map_rep1_S81, run3/CPA_Map_rep3_S80) and the real hg38 TTAA
# reference: with V3 (this line, tt=Start in stage 04), minus-strand hit rate is ~61-69%,
# matching the plus strand's own ~58-70% baseline; with V3-1 and the "matching" tt=Start+1,
# minus-strand hit rate collapses to ~1.6-1.9%. This also matches both igoru's own current
# production script (lidiag/R/scripts/analyzeTRIPMapping_uniform.R) and the original,
# unmodified ArtyomTyukaev/trip-pipeline script -- both use plain V3, never V3-1.
tabStarts$Start = ifelse(tabStarts$V6=="+",tabStarts$V2,tabStarts$V3)
tabStarts$Length = tabStarts$V3-tabStarts$V2
tabBarcodes$V1=substr(tabBarcodes$V1,2,1000)
print("Merging by read name")
tabStarts$Name = sapply(strsplit(as.character(tabStarts$V4), "/"), "[", 1)
length(intersect(tabStarts$Name,tabBarcodes$V1))
merged=merge(tabBarcodes,tabStarts,by.x=1,by.y="Name")
#table(nchar(merged$V3.x))

# Keep only barcodes of the expected length
merged = subset(merged,nchar(merged$V3.x)==barcodeLen)

# Remove unneeded columns
merged$V1=NULL
merged$V4=NULL
merged$V9=NULL
dim(merged)

nBefore = nrow(merged)
# Remove redundant
merged=unique(merged)
nAfter = nrow(merged)
print(paste("Removing redundant rows reduced from ",nBefore,"to",nAfter))



# Remove plasmids with globin integration sequence
isGlobin = merged$V1.y=="chr11" & merged$Start>=5220000 & merged$Start<=5240000
print(paste(sum(isGlobin),"Globin integrations will be excluded"))
merged=subset(merged,!isGlobin)
#sort(table(paste0(merged$V1.y,":",merged$Start)),decreasing = T)[1:10]

print("Aggregating")

spChr = split(merged,merged$V1.y)

agg=list()
for (chr in names(spChr))
{
  if (sum(grepl("alt|random",chr))>0)
    next;
  agg[[chr]] = aggregate(spChr[[chr]][,c("V2.x","V3.x")],by=spChr[[chr]][,c("V1.y","Start","V6","Length")],function(x){paste0(x,collapse=",")})
}

head(agg$chr2)
agg$chr2[agg$chr2$Start==32916313,]

screenCount <- function(x){
  umis=unlist(strsplit(x[1],","));
  barcodes=unlist(strsplit(x[2],","));
  if (length(umis)!=length(barcodes))
  {
    print("Discordant!")
    print(x[3])
    print(x[4])
    print(x[5])
    print(umis)
    print(barcodes)
  }
  sp=split(umis,barcodes);
  counts=sapply(sp,length);
  myMax = max(counts)
  sum(counts>=myMax*0.01)
}

screenFunc <- function(x){
                          umis=unlist(strsplit(x[1],","));
                          barcodes=unlist(strsplit(x[2],","));
                          sp=split(umis,barcodes);
                          counts=sapply(sp,length);
                          sp=sp[order(counts,decreasing = T)]
                          counts=counts[order(counts,decreasing = T)]
                          outVec = list()
                          outVec[[names(sp)[1]]]=length(sp[[1]])
                          maxCount = max(sapply(sp,length))
                          if (sum(counts>=maxCount*0.01)>=50)
                          {
                            print(paste("Noisy position",x[3],x[4]))
                            ""
                          } else {
                            if (length(sp)>1)
                            {
                              for (i in 2:length(sp))
                              {
                                myBarcode = names(sp)[i]
                                myCount = length(sp[[i]])
                                if (myCount < 0.01*maxCount)
                                {
                                  next;
                                }
                                for (j in 1:length(outVec))
                                {
                                  if (StrDist(myBarcode,names(outVec)[j])<=3)
                                  {
                                    outVec[[j]]=outVec[[j]]+myCount
                                  } else {
                                    outVec[[names(sp)[i]]]=myCount
                                  }
                                }
                              }
                            }
                            paste0(paste0(names(outVec),":",unlist(outVec)),collapse=",")
                          }
                        }
#agg$umiNames="Globin"
print("Finished aggregating")
print("Starting to filter")
#notGlobin = agg$V1.y!="chr6" | agg$Start<5220000 | agg$Start>5240000
for (chr in names(agg))
{
  print(chr)
  agg[[chr]]$umiNames=apply(agg[[chr]][,c("V2.x","V3.x","V1.y","Start")],1,screenFunc)  
  agg[[chr]]=subset(agg[[chr]],agg[[chr]]$umiNames!="")
}

table(agg$chr2$Length)

aggAll = do.call("rbind",agg)
dim(aggAll)
agg=aggAll


print("Finished screenFunc")

print("Buidling output")
df= data.frame(Chrom=agg$V1.y,Start=agg$Start,End=agg$Start+1,Name=agg$umiNames,Score=1000,Strand=agg$V6,Length=agg$Length)
df$counts = unlist(lapply(strsplit(df$Name,","),length))
df$cleanCounts = 0

df$Barcodes=unlist(lapply(lapply(strsplit(df$Name,","),substr,1,barcodeLen),paste0,collapse=","))

#sapply(strsplit(df[df$counts<5,"Name"],","),Rncc::RcppNumberConnectedComponents,2)
#df$cleanCounts[df$counts<100] = unlist(lapply(strsplit(df[df$counts<100,"Name"],","),function(x){Rncc::RcppNumberConnectedComponents(x,2)}))
df$cleanCounts = sapply(df$counts,min,10000)
write.table(df,paste0(path,name,".barcodes_poses.txt"),col.names=F,row.names=F,quote=F,sep="\t")
df$Name = substr(df$Name,1,200)
print("Writing")
dfBED = df[,c("Chrom","Start","End","Name","cleanCounts","Strand")]
options(scipen=10)
write.table(dfBED[order(dfBED$Chrom,dfBED$Start),],paste0(path,name,".barcodes.bed"),col.names=F,row.names=F,quote=F,sep="\t")
