  library(doMC)
  library(ggplot2)
  registerDoMC(cores = 4)
  library(data.table)
  library(dplyr)
  library(plyr)
  library(stringr)
  library(foreach)
  library(doParallel)
  #set number of cores
  registerDoParallel(cores =4)
#Read file
yyy = fread(args[1],header=FALSE,stringsAsFactors = FALSE,quote = "")
#colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')
rmap <- fread(args[2],header=FALSE,stringsAsFactors = FALSE,quote = "")
colnames(rmap) <- c('chrom','start','end','Frag_id')
baitmap <- fread(args[3],header=FALSE,stringsAsFactors = FALSE,quote = "")
colnames(baitmap) <- c('chrom','start','end','Frag_id','bait')
rmap <- rmap[rmap$Frag_id < 143731]
rsplit <- split(rmap,rmap$chrom)




#Read size file
chrsize <- fread("~/Project/Data/NiGa/D13.sizes",header=FALSE,stringsAsFactors = FALSE,quote = "")
chrsize <- chrsize[order(as.numeric(as.character(chrsize$V1)))]

##############
chkl <- split(yyy,yyy$dnachrom)
rm(yyy)
gc()
bin = NULL
for (i in 1:22){
  bin[[i]] <- transform(chkl[[i]], group = cut(dnachromStart,
                                           breaks=seq(from = 0, to = chrsize$V2[i], by = 20000 )))
}
rm(chkl)
gc(reset = TRUE)



for (i in 1:length(bin)) {
  bin[[i]] <-transform(bin[[i]], Freq = ave(seq(nrow(bin[[i]])),group, FUN = length))
  bin[[i]] <- bin[[i]][,-c(2:3)]
}

#Defining function to check integer(0)
is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

#getting otherendID column and distal length

k = NULL


foreach (k = 1:22) %:% 
  foreach (i = 1:length(bin[[k]]$rnachrom)) %dopar%  {
    bin[[k]]$otherID[i] <- ifelse(is.integer0(rsplit[[k]]$Frag_id[which(rsplit[[k]]$start <= bin[[k]]$dnachromStart[i] & rsplit[[k]]$end >= bin[[k]]$dnachromStart[i])]),NA,rsplit[[k]]$Frag_id[which(rsplit[[k]]$start <= bin[[k]]$dnachromStart[i] & rsplit[[k]]$end >= bin[[k]]$dnachromStart[i])])
    bin[[k]]$distLen[i] <- ifelse(bin[[k]]$rnachrom != k, NA,min(baitmap$start[which(bin[[1]]$rnaID[i] == baitmap$bait)]-bin[[1]]$dnachromStart[i],baitmap$start[which(bin[[1]]$rnaID[i] == baitmap$bait)]-bin[[1]]$dnachromEnd[i],baitmap$end[which(bin[[1]]$rnaID[i] == baitmap$bait)]-bin[[1]]$dnachromStart[i],baitmap$end[which(bin[[1]]$rnaID[i] == baitmap$bait)]-bin[[1]]$dnachromEnd[i]) )	
    
}
#Defining otherendLength column
for (i in 1:length(rmap$chrom)){
  if (!(rmap$Frag_id[i] %in% baitmap$Frag_id)){
    rmap$otherLen[i] <- rmap$end[i]-rmap$start[i]
  } else {
    rmap$otherLen[i] <- NA
  }
}

all_chr_freq <- bind_rows(bin, .id = "column_label")
all_chr_freq <- all_chr_freq[with(all_chr_freq,order(as.numeric(column_label),as.numeric(group))),]


all_chr_freq$baitID <- baitmap$Frag_id[match(all_chr_freq$rnaID,baitmap$bait)]
all_chr_freq$otherLen <- rmap$otherLen[match(all_chr_freq$otherID,rmap$Frag_id)]







chinput <- data.frame(all_chr_freq$baitID,all_chr_freq$otherID,all_chr_freq$Freq,all_chr_freq$otherLen,all_chr_freq$distLen)
chinput <- chinput[with(chinput,order(as.numeric(chinput$all_chr_freq.baitID))),]
chinput <- chinput[!is.na(chinput$all_chr_freq.otherID),]
chinput <- chinput[!is.na(chinput$all_chr_freq.baitID),]
chinput[sapply(chinput, is.infinite )] <- NA
chinput$all_chr_freq.distLen <- as.numeric(chinput$all_chr_freq.distLen)
#Define the output file
write.table(chinput,"~/Project/Data/1.chinput",sep = "\t",col.names = FALSE,row.names = FALSE)
