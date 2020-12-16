library(doMC)
library(ggplot2)
registerDoMC(cores = 4)
library(data.table)
library(dplyr)
library(plyr)
library(stringr)
#Read file
yyy = fread("~/Project/Data/NiGa/Day13_done.bedpe",header=FALSE,stringsAsFactors = FALSE,quote = "")
#colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')





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



chr_freq = NULL
test = NULL
t = NULL
for (i in 1:length(bin)){
  #bin[[i]]$group <- paste(bin[[i]]$rnachrom,bin[[i]]$group)
  #test[[i]] <- data.frame(bin[[i]]$rnaID,bin[[i]]$group)
  bin[[i]] <-transform(bin[[i]], Freq = ave(seq(nrow(bin[[i]])),group, FUN = length))
  #t[[i]] <- as.data.frame.matrix(table(test[[i]]))
}

for  (i in 1:22){
  bin[[i]] <- bin[[i]][,-c(2:3)]
}


all_chr_freq <- bind_rows(chr_freq, .id = "column_label")
all_chr_freq <- arrange(all_chr_freq, group)
all_chr_freq <- all_chr_freq[all_chr_freq$rnaID != "Intergenic"]
#t <- as.data.frame.matrix(table(test))
all_chr_freq <- all_chr_freq[with(all_chr_freq,order(as.numeric(column_label),as.numeric(group))),]





##########
rmap <- fread("~/Project/Data/Des/D13.rmap",header=FALSE,stringsAsFactors = FALSE,quote = "")
colnames(rmap) <- c('chrom','start','end','Frag_id')
baitmap <- fread("~/Project/Data/Des/D13.baitmap",header=FALSE,stringsAsFactors = FALSE,quote = "")
colnames(baitmap) <- c('chrom','start','end','Frag_id','bait')

#Defining function to check integer(0)
is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

#getting otherendID column
rsplit <- split(rmap,rmap$chrom)
otherID = NULL
for (k in 1:22){
  for (i in 1:length(bin[[k]]$rnachrom)){
    if (rsplit[[k]]$Frag_id < 143730){
      print(i)
      if (is.integer0(rsplit[[k]]$Frag_id[which(rsplit[[k]]$start <= bin[[k]]$dnachromStart[i] & rsplit[[k]]$end >= bin[[k]]$dnachromStart[i])])){
        bin[[k]]$otherID[i] <- NA
    } else{
        bin[[k]]$otherID[i] <- rsplit[[k]]$Frag_id[which(rsplit[[k]]$start <= bin[[k]]$dnachromStart[i] & rsplit[[k]]$end >= bin[[k]]$dnachromStart[i])]
    }
  }
}




all_chr_freq$intg <- as.character(all_chr_freq$group)
all_chr_freq$intg = gsub("\\]","",all_chr_freq$intg)
all_chr_freq$intg = gsub("\\(","",all_chr_freq$intg)

#a = str_split_fixed(all_chr_freq$intg,",",2)
foo <- data.frame(do.call('rbind',strsplit(all_chr_freq$intg,',',fixed = TRUE)))
ID = NULL
N = NULL
foo$X3 <- all_chr_freq$column_label
foo$X1 <- as.character(foo$X1)
foo$X2 <- as.character(foo$X2)
foo$X1 <- as.numeric(foo$X1)
foo$X2 <- as.numeric(foo$X2)
#for (j in 1:length(rmap$chrom)){
#  for (i in 1:length(all_chr_freq$column_label)){
#    if ((rmap$chrom[j] == foo$X3[i]) && (between(rmap$start[j],foo$X1[i],foo$X2[i])) && (!(rmap$Frag_id[j] %in% baitmap$Frag_id))){
#      print(rmap$Frag_id[j])
#      all_chr_freq$otherID[i] <-rmap$Frag_id[j]
#    }
#  }
#}
all_chr_freq$baitID <- baitmap$Frag_id[match(all_chr_freq$rnaID,baitmap$bait)]
all_chr_freq$otherLen <- rmap$otherLen[match(all_chr_freq$otherID,rmap$Frag_id)]


baitmap$bait[1]




for (i in 1:length(foo$X1)){
  for (j in 1:22){
    if (rmap$chrom == j) {
      which(rmap$start >= foo$X1[i] & rmap$start <= foo$X2[i])
    }
  }
}
  
#turn factor to character , uuse genomic ranges. 

#Defining otherendLength column
for (i in 1:length(rmap$chrom)){
  if (!(rmap$Frag_id[i] %in% baitmap$Frag_id)){
    rmap$otherLen[i] <- rmap$end[i]-rmap$start[i]
  } else {
    rmap$otherLen[i] <- NA
  }
}

#Distal column
for (k in 1:22){
  for (i in 1:length(bin[[k]]$rnachrom)){
    print(i)
    if (bin[[k]]$rnachrom[i] != bin[[k]]$dnachrom[i]){
      bin[[k]]$distLen[i] <- NA
  } else {
        bin[[k]]$distLen[i] <- min(baitmap$start[which(bin[[k]]$rnaID[i] == baitmap$bait)]-bin[[k]]$dnachromStart[i],baitmap$start[which(bin[[k]]$rnaID[i] == baitmap$bait)]-bin[[k]]$dnachromEnd[i],baitmap$end[which(bin[[k]]$rnaID[i] == baitmap$bait)]-bin[[k]]$dnachromStart[i],baitmap$end[which(bin[[k]]$rnaID[i] == baitmap$bait)]-bin[[k]]$dnachromEnd[i])
    }
  }  
}

write.table(all_chr_freq,"~/Project/Data/D13.chinput",sep = "\t")

chinput <- data.frame(all_chr_freq$baitID,all_chr_freq$otherID,all_chr_freq$Freq,all_chr_freq$otherLen,all_chr_freq$distLen)
chinput <- chinput[with(chinput,order(as.numeric(chinput$all_chr_freq.baitID))),]
chinput <- chinput[!is.na(chinput$all_chr_freq.otherID),]
chinput <- chinput[!is.na(chinput$all_chr_freq.baitID),]
chinput[sapply(chinput, is.infinite )] <- NA
chinput$all_chr_freq.distLen <- as.numeric(chinput$all_chr_freq.distLen)
chinput
write.table(chinput,"~/Project/Data/1.chinput",sep = "\t")