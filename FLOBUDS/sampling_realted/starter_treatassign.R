rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
setwd("~/Documents/git/proterant/fLOBUDS")


d<-read.csv("flobud_datapoints.csv") ###read in data

d$classification<-paste(d$Genus,d$species,sep="_") ###make a column contains genus_species

###how many samples need from each individual?
table(d$classification)
#8 treatments, 6 cuttings per treatment within species= 48 total cuttings per pecies
cuttingper<-as.data.frame(table(d$classification))
cuttingper$cutting_number<-48/cuttingper$Freq

cuttingper<-filter(cuttingper, Freq>6) ### This sheets says how many cuttings per

#assign treatment
#sample(seq(from = 1, to = 8, by = 1), size = 8, replace = FALSE)
#adapted from Dan F's code
chill <- gl(2, 1, labels = c("N","Y"))
force <- gl(2, 2, length=4, labels = c("C","W"))
photo<-gl(2,4,label=c("S","L"))
treatcode <- paste(substr(chill, 1, 4), substr(force, 1, 4), substr(photo, 1, 4),sep="_") 

treat <- data.frame(chill,force,photo,treatcode)

###This is where I leave off, can't make sense of the loop

dx <- vector()

for(i in 1:nrow(d) ) { # i = 1
  
  
  xx <- paste(d[i,"Individual"], formatC(1:8, width = 2, flag = "0"), sep = "_")
  
  xx <- data.frame(xx)
  
  xx$sp <- substr(xx[,1], 1, 6)
  xx$rep <- substr(xx[,1], 7, 8)
  xx$ind <- substr(xx[,1], 1, 8)
  
  
  names(xx)[1] = "id"
  
  # Assign treatments. Randomize rows of treatment dataframes and apply to this individual
  xx <- data.frame(xx, treat[sample(1:nrow(treat)),])
  
  dx <- rbind(dx, xx)
}


###try dplyer
d<-d %>% group_by(classification)