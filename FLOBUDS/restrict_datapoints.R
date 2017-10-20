###THis file updates the data point to only include species I want

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

setwd("~/Documents/git/proterant/FLOBUDS")
library(dplyr)

dat<-read.csv("flobud_datapoints.csv",header=TRUE) 
dat$nomen<-paste(dat$Genus,dat$species,sep="_")
unique(dat$nomen)
do<-filter(dat,nomen!="Viburnum_cassinoides")
do<-filter(do,nomen!="Fagus_grandifolia")
do<-filter(do,nomen!="Quercus_alba")
do<-filter(do,nomen!="Quercus_rubra")
do<-filter(do,nomen!="Fraxinus_americana")
do<-filter(do,nomen!="Lindera_benzoin")


remove<-c("Viburnum_cassinoides","Fagus_grandifolia","Quercus_alba","Quercus_rubra","Fraxinus_americana","Lindera_benzoin")
unique(do$nomen)
      
write.csv(do,"final_flobud_datapoints.csv",row.names = FALSE)
