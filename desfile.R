library(doMC)
library(ggplot2)
registerDoMC(cores = 4)
library(data.table)
library(dplyr)
library(plyr)
library(compare)

yyy = fread("~/Project/Data/NiGa/Day13_with.bedpe",header = FALSE,stringsAsFactors = FALSE,quote = "")
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')
#Make length colums
yyy$len <- abs((yyy$dnachromStart+yyy$dnachromEnd)/2-(yyy$rnachromStart+yyy$rnachromEnd)/2)

chrsize <- fread("~/Project/Data/NiGa/D13.sizes",header=FALSE,stringsAsFactors = FALSE,quote = "")


# Reading chr.sizes
chrsize <- chrsize[order(as.numeric(as.character(chrsize$V1)))]


############ Generating rmap file
lit1 = NULL
lit2 = NULL
lit3 = NULL
rrmap = NULL
for (i in 1:22){
  lit1[[i]] <- seq(1,chrsize$V2[i]-20000, by = 20001)
  lit2[[i]] <- seq(20001,chrsize$V2[i], by = 20001)
  lit3[[i]] <- rep(i,length(lit1[[i]]))
  rrmap[[i]] = data.frame("chr" = lit3[[i]],"start" = lit1[[i]], "end" = lit2[[i]])
}



df <- bind_rows(rrmap, .id = "column_label")
df$column_label <- NULL



df

###############3

######
#Reading the transript file
unqID = unique(yyy$rnaID)
IDs = fread("~/Project/Code/trans.bedpe",header=FALSE,stringsAsFactors = FALSE,quote = "")
#### Generating trial baitmap file


data = filter(IDs,V4 %in% unqID)
data = data %>% mutate_all(~gsub("chr","",.))

data[,sapply(V1,is.numeric)]
data = data[!data$V1 %in% c("X","Y"),]

data <- data[with(data,order(as.numeric(V1,V2))),]


#Naming columns in baitmap
colnames(data) <- c('chr','start','end','baitAnnotation')
######

bait_try <- data.frame(data[,1:3])

rmap_try <- rbind(df,bait_try)
rmap_try <- rmap_try[with(rmap_try,order(as.numeric(chr),as.numeric(start))),]

rmap_try[,c(1,2,3)] <- sapply(rmap_try[,c(1,2,3)], as.numeric)
k = (length(rmap_try$chr)-2)
for (i in 1:k){
  if ((rmap_try$end[i]-rmap_try$start[i]==20000)&&(rmap_try$end[i] > rmap_try$start[i+1])){
    rmap_try[i,3] <- rmap_try[i+1,2]-1
    #rmap_try$end[i] <- rmap_try$start[i+1]-1
    print(i)
  }
}
rmap_try
rmap_try$fragmentID <- 1:nrow(rmap_try)
jis = 0
data$fragmentID <- NULL
data[,c(1,2,3,4)] <- sapply(data[,c(1,2,3,4)], as.numeric)
data <- data[with(data,order(as.numeric(chr),as.numeric(start))),]
for (i in 1:length(rmap_try$start)){
  for (j in 1:length(data$start)){
    if ((rmap_try$end[i]-rmap_try$start[i] != 20000) &&(rmap_try$start[i] == data$start[j])){
      data$fragmenID[j] = rmap_try$fragmentID[i]
      #jis = jis+1
    }
  }
}

write.table(rmap_try,"~/Project/Data/D13.rmap",sep = "\t")



file1 = fread("~/Project/Data/D13.rmap",header = FALSE,stringsAsFactors = FALSE,quote = "")
colnames(file1) <- c('number','chr','start','end','fragmentID')
file1 <- as.data.frame(sapply(file1, function(number) gsub("\"","",number)))
file1$fragmentID <- file1$number
file1$number <- NULL
file1$fragmentID <- as.numeric(as.character(file1$fragmentID))
file1$chr <- as.numeric(as.character(file1$chr))
file1$start <- as.numeric(as.character(file1$start))
file1$end <- as.numeric(as.character(file1$end))
fragID = NULL
for (i in 1:length(file1$chr)){
  if (file1$fragmentID[i] > length(df$chr)){
    fragID[i] = file1$fragmentID[i]
  }
}
fragID <- fragID[!is.na(fragID)]
data$fragmentID <- fragID
file1 <- file1[(file1$end != 0),]
file1 <- file1[file1$fragmentID != 1e+05,]
for (i in 1:length(file1$chr)){
  if ((abs(file1$start[i+1]-file1$end[i]) > 1) && (file1$fragmentID[i+1] < 143730)){
    print(i)
    file1$start[i+1] = file1$end[i]+1
  }
}
file1 <- file1[!(file1$start > file1$end),]

data <- data[,c(1,2,3,5,4)]

data <- as.data.frame(sapply(data, function(baitAnnotation) gsub("\"","",baitAnnotation)))

data
write.table(data,"~/Project/Data/D13.baitmap",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(file1,"~/Project/Data/D13.rmap",sep = "\t",col.names = FALSE,row.names = FALSE)

#file1[8:15,1:5]



#for (i in 1:length(rmap_try$start)){
 # if ( file1$number[i] > 143731 ){
  #  file1$fragmentID[i] = file1$number[i]
  #}
#}

#try = file1 %>% filter(number > 143731)
#for (i in 1:length(try$number)){
 # try$baitAnnotation[i] = data$baitAnnotation[i]
#}
#try
################

windowSize = 20000
map = NULL
for (whichChr in unique(yyy$rnachrom)){
  start = seq(from=0,to=chrsize[chrsize$V1==whichChr]$V2,by=windowSize)
  end = start+windowSize
  some[whichChr] = rep(whichChr,length(start))
  map[[whichChr]] = data.frame(some[whichChr],start,end)
}
map[[2]]
for (i in 1:22){
  map[[i]] <- transform(map[[i]], group = cut(end,
                                 breaks=seq(from = 0, to = chrsize$V2[i], by = 20000 )))
}
df <- bind_rows(map, .id = "column_label")
df$column_label <- NULL
colnames(df)[1] <- 'chr'


#write.table(data,"~/Project/Data/D13.baitmap",sep = "\t")
#write.table(rmap_try,"~/Project/Data/D13.rmap",sep = "\t")


