bed <- read.csv("~/Project/Data/NiGa/spliaa.bedpe",header=FALSE,sep="\t",stringsAsFactors=FALSE, quote="")
library(doMC)
library(ggplot2)
registerDoMC(cores = 4)
library(data.table)
library(dplyr)
library(plyr)
yyy = fread("~/Project/Data/NiGa/uniq.TGACCA_D13_3.bedpe",header=FALSE,sep = "\t",stringsAsFactors = FALSE,quote = "")
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
yyy$len <- abs((yyy$dnachromStart+yyy$dnachromEnd)/2-(yyy$rnachromStart+yyy$rnachromEnd)/2)
mean(yyy$len)
yyy$match <- NULL
yyy[yyy$rnachrom == yyy$dnachrom,"match"] <-"1"
yyy[yyy$rnachrom != yyy$dnachrom,"match"] <-"0"
vc <- '1'
schr <- yyy[is.element(yyy$match,vc),]
table(schr$len)
barplot(table(schr$len))
ggplot(schr,aes(x = len))+
  geom_histogram(aes(y = stat(count)/sum(count)))+
  scale_y_continuous(labels = scales::percent)

lsd <- c(unique(schr$rnaID))
lsd[2]

schr[order(schr$rnaID,decreasing = TRUE)]
sdf <- daply(schr,)
kls <- split(schr,schr$rnaID)
gir <- NULL
for (i in 1:length(kls)){
  gir[i] <- length(kls[[i]]$rnachrom)
}
#########################
pdf(file="Plots.pdf")
for (i in 1:length(kls)){
  if (length(kls[[i]]$rnachrom)>100){
    ggplot(kls[[i]],aes(x = len))+
      geom_histogram(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)
  }
}
dev.off()
###################
plot_list = list()
for (i in 1:length(kls)){
  if (length(kls[[i]]$rnachrom)>100){
    p = ggplot(kls[[i]],aes(x = len))+
      geom_histogram(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)
  } else {
    p = ggplot(kls[[i]],aes(x = len))+
      geom_freqpoly(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)
  }
  plot_list[[i]] = p
}
tt <- NULL
for (i in 1:length(gir)){
  if (gir[i]>100){
    tt[i] <- 1
  }
}
count(tt)
class(gir)
for (i in 1:length(kls)){
  file_name = paste("iris_plot_", i, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}
pdf("plots.pdf")
for (i in 1:length(kls)){
  print(plot_list[[i]])
}
dev.off
print(plot_list[[3]])

#######


ggplot(kls[[1]],aes(x = len))+
  geom_freqpoly(aes(y = stat(count)/sum(count)))+
  scale_y_continuous(labels = scales::percent)

plot(kls[[2]]$len)
frequency()
for (i in 1:length(lsd)){
  f_df[[i]] <- filter(schr,rnaID == lsd[i] )
  ggplot(f_df[[i]],aes(x = len))+
    geom_histogram(aes(y = stat(count)/sum(count)))+
    scale_y_continuous(labels = scales::percent)
}
#---------------------------------------------#

plot_list = list()
for (i in 1:length(kls)){
  if (length(kls[[i]]$rnachrom)>100){
    p = ggplot(kls[[i]],aes(x = len))+
      geom_histogram(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)
  }
  plot_list[[i]] = p
}
for (i in 1:3051){
  file_name = paste("iris_plot_", i, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}
pdf("plots.pdf")
for (i in 1:3051){
  print(plot_list[[i]])
}
dev.off
