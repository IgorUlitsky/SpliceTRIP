library(DescTools)
library(dplyr)
library(purrr)
library(tidyr)

path="/home/labs/ulitsky/shared/data/human/TRIP/set1_combined_1/" 
nameMap = "CPA_noP.txt"
barcodeLen = 18

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path = myArgs[1]
  nameMap = myArgs[2]
  # myArgs[3] is optional (barcodeLen keeps its default of 18 above when omitted) -- indexing
  # a shorter myArgs vector out of bounds returns NA rather than erroring, and as.integer(NA)
  # is silently NA, not 18. That NA then propagates through nchar(barcodes)==barcodeLen
  # below as an all-NA comparison, which silently empties out every position's barcode list
  # with no error at all (confirmed in the sibling get_genuine_barcodes_IntronTRIP.R script:
  # running with only path+nameMap produced a "successful" run with fully empty output).
  if (length(myArgs)>=3) barcodeLen = as.integer(myArgs[3])
}
print(path)
print(nameMap)

extract_info <- function(barcode_counts) {
  # Split the string into barcodes and their counts
  parts <- strsplit(barcode_counts, ",")[[1]]
  
  if(length(parts) == 0) return(list(max_barcode = NA, max_count = NA, remaining_barcodes = NA))  # Handle empty input
  
  df_barcodes <- do.call(rbind, strsplit(parts, ":"))
  colnames(df_barcodes) <- c("barcode", "count")
  df_barcodes <- as.data.frame(df_barcodes, stringsAsFactors = FALSE)
  df_barcodes$count <- as.numeric(df_barcodes$count)
  
  # Aggregate counts for each barcode
  df_barcodes <- df_barcodes %>%
    group_by(barcode) %>%
    summarise(count = sum(count), .groups = 'drop')
  
  # Find the barcode with the maximum count
  max_row <- df_barcodes[which.max(df_barcodes$count), ]
  max_barcode <- max_row$barcode
  max_count <- max_row$count
  
  # Find the remaining barcodes and their counts
  remaining_barcodes <- df_barcodes %>%
    filter(barcode != max_barcode) %>%
    mutate(remaining = paste(barcode, count, sep = ":")) %>%
    pull(remaining) %>%
    paste(collapse = ",")
  
  if (remaining_barcodes == "") remaining_barcodes <- NA
  
  return(list(max_barcode = max_barcode, max_count = max_count, remaining_barcodes = remaining_barcodes))
}

# Function to process the dataset
process_dataset <- function(data) {
  data <- data %>%
    rowwise() %>%
    mutate(
      info = list(extract_info(barcode))
    ) %>%
    ungroup() %>%
    mutate(
      actual_barcode = sapply(info, "[[", "max_barcode"),
      actual_count = sapply(info, "[[", "max_count"),
      remaining_barcodes = sapply(info, "[[", "remaining_barcodes")
    ) %>%
    select(-info)  # Remove the temporary 'info' column
  
  return(data)
}

# Function to calculate Hamming distance
hamming_distance <- function(str1, str2) {
  stringdist::stringdist(str1, str2, method = "hamming")
}

# Function to process the dataset and update actual_count based on Hamming distance
update_actual_count <- function(data) {
  # Iterate over each row
  for (i in 1:nrow(data)) {
    if (!is.na(data$remaining_barcodes[i])) {
      # Split remaining_barcodes into individual barcodes and counts
      remaining_barcodes_list <- strsplit(data$remaining_barcodes[i], ",")[[1]]
      additional_count <- 0
      
      # Iterate over each remaining barcode
      for (barcode_count in remaining_barcodes_list) {
        bc <- strsplit(barcode_count, ":")[[1]]
        barcode <- bc[1]
        count <- as.numeric(bc[2])
        
        # Calculate Hamming distance and update additional_count if distance <= 2
        if (hamming_distance(data$actual_barcode[i], barcode) <= 2) {
          additional_count <- additional_count + count
        }
      }
      
      # Update actual_count
      data$actual_count[i] <- data$actual_count[i] + additional_count
    }
  }
  
  return(data)
}
process_other_barcodes <- function(df) {
  
  df <- df %>% arrange(desc(actual_count))
  
  # Vector for storing barcodes that need to be deleted
  to_remove <- vector("logical", nrow(df))
  counter <- 0

  for (i in 2:nrow(df)) {
    for (j in 1:(i - 1)) {
      # Hamming distance check
      if (hamming_distance(df$actual_barcode[i], df$actual_barcode[j]) <= 2) {
        # Marking a row for deletion
        to_remove[i] <- TRUE
        break
      }
    }
    counter <- counter + 1
    if (counter %% 1000 == 0) {
      print(paste(counter, ' rows done'))
    }
  }
  
  # Return filtered dataset
  df <- df[!to_remove, ]
  
  return(df)
}

# Read the data
print('Reading data')


datas = list()
for (fileName in read.table(paste0(path, nameMap))$V1)
{
  if (file.size(paste0(path,fileName))>0)
  {
    data <- read.table(paste0(path,fileName), sep = "\t", stringsAsFactors = FALSE)
    # Set column names
    colnames(data) <- c('chr', 'start', 'end', 'barcode', 'count', 'strand')
    datas[[fileName]]=data
  }
}


bound = do.call(rbind,datas)
dim(bound)
aggregated = aggregate(bound$barcode,by=bound[,c(1:3,6)],FUN=paste0,collapse=",")
aggregated$count=sapply(aggregated$x,function(a){length(strsplit(a,",")[[1]])})

aggregated=aggregated[,c("chr","start","end","x","count","strand")]
names(aggregated) <- c('chr', 'start', 'end', 'barcode', 'count', 'strand')
dim(aggregated)

filterBarcodes<-function(x)
{
  chr=x[1]
  start=x[2]
  end=x[3]
  strand=x[6]
  myStr=x[4]
  pairs = sapply(strsplit(myStr,","),strsplit,":")
  pairs=pairs[sapply(pairs,length)==2]
  barcodes = sapply(pairs,"[",1)
  counts = sapply(pairs,"[",2)
  counts=counts[nchar(barcodes)==barcodeLen & !grepl("N",barcodes)]
  barcodes=barcodes[nchar(barcodes)==barcodeLen & !grepl("N",barcodes)]
  countMat = data.frame(barcode=barcodes,count=counts)
  if (nrow(countMat)>0)
  {
    countMat = aggregate(as.numeric(countMat$count),by=list(countMat$barcode),FUN=sum)
    countMat=countMat[order(countMat$x,decreasing = T),]
    recontedMat = countMat
    if (nrow(countMat)>1)
    {
      keep = rep(T,nrow(recontedMat))
      for (i in 2:nrow(recontedMat))
      {
        for (j in 1:(i-1))
        {
          if (StrDist(recontedMat[i,"Group.1"],
                      recontedMat[j,"Group.1"])<=2)
          {
            recontedMat[j,"x"]=recontedMat[i,"x"]+recontedMat[j,"x"]
            keep[i]=F
            break;
          }
        }
      }
      recontedMat=recontedMat[keep,]
    }
    data.frame(cbind(chr,start,end,recontedMat,strand))
  }
}


x=unlist(aggregated[1,])
filteredBarcodes = apply(aggregated,1,filterBarcodes)
filtered=do.call("rbind",filteredBarcodes)
dim(filtered)
locationCounts = table(filtered$Group.1)

singletons = subset(filtered,locationCounts[filtered$Group.1]==1)
dim(singletons)
write.table(singletons[order(singletons$chr,as.numeric(singletons$start)),], paste0(path,nameMap,".retained_barcodes_genuine.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
forUCSC = singletons[order(singletons$chr,as.numeric(singletons$start)),]
forUCSC[,5]=sapply(forUCSC[,5],min,1000)
write.table(forUCSC, paste0(path,nameMap,".retained_barcodes_genuine_forUCSC.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

multiMappers = subset(filtered,locationCounts[filtered$Group.1]>1)
dim(multiMappers)
write.table(multiMappers[order(multiMappers$chr,as.numeric(multiMappers$start)),], paste0(path,nameMap,".retained_barcodes_multimappers.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
forUCSC = multiMappers[order(multiMappers$chr,as.numeric(multiMappers$start)),]
forUCSC[,5]=sapply(forUCSC[,5],min,1000)
write.table(forUCSC, paste0(path,nameMap,".retained_barcodes_multimappers_forUCSC.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

### Lydia's old code

#print('Selecting the barcode with the maximum score')
# Selects the barcode with the maximum score
#processed_data <- process_dataset(aggregated)

#print("Increasing a barcode's score if other barcodes mapped to that location are its mutants")
# Increases a barcode's score if other barcodes mapped to that location are its mutants 
#processed_data_new <- update_actual_count(processed_data)
#dim(processed_data_new)

#processed_data_new$remaining_barcodes

#print('Removing barcodes mutants with minimal score at other positions in the genome')
# Removes barcodes mutants with minimal score at other positions in the genome
#processed_data_new_new <- process_other_barcodes(processed_data_new)
#processed_data_new_new=processed_data_new

print('Saving results')
# Saving results
#write.table(processed_data_new_new, paste0(path,nameMap,".retained_barcodes_genuine.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
