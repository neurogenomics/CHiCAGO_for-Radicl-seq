#import libraries
library(doMC)
library(ggplot2)
registerDoMC(cores = 4)
library(data.table)
library(dplyr)
library(plyr)
#Read file
yyy = fread("~/Project/Data/NiGa/uniq.TGACCA_D13_3.bedpe",header=FALSE,sep = "\t",stringsAsFactors = FALSE,quote = "")
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')

#Make length colums
yyy$len <- abs((yyy$dnachromStart+yyy$dnachromEnd)/2-(yyy$rnachromStart+yyy$rnachromEnd)/2)
mean(yyy$len)

#Remove trans inetractions
yyy$match <- NULL
yyy[yyy$rnachrom == yyy$dnachrom,"match"] <-"1"
yyy[yyy$rnachrom != yyy$dnachrom,"match"] <-"0"
vc <- '1'
schr <- yyy[is.element(yyy$match,vc),]
#########

#Sort the file
kls <- split(schr,schr$rnaID)

######################### Get the plots
###################
setwd("~/Project/Plots/D13_3/")

plot_list = list()
for (i in 1:length(kls)){
  if (length(kls[[i]]$rnachrom)>100){
    p = ggplot(kls[[i]],aes(x = len))+
      geom_histogram(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)+
      ggtitle(kls[[i]]$rnaID[1])
    plot_list[[i]] = p
  } else {
    p = ggplot(kls[[i]],aes(x = len))+
      geom_freqpoly(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)+
      ggtitle(kls[[i]]$rnaID[1])
    plot_list[[i]] = p
  }
}
for (i in 1:length(kls)){
  file_name = paste("distance_plot_", kls[[i]]$rnaID[1], ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}
pdf("plots.pdf")
for (i in 1:length(kls)){
  print(plot_list[[i]])
}
dev.off


#######



