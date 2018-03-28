##basic exploration of harvard forest data
###offset (time between 1 phenophase and the other) would give you the ability to predict loss of hysteranthy once sensativities have been calculated
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(plotrix)
library(gdata)
library(nlme)
library(scales)
library(arm)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/git/proterant/Data")
hf<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
mean_hf<-read.csv("hf003-06-mean-spp.csv",header = TRUE)
hf<-mutate(hf, buds_seq = bb.jd -fbb.jd)
hf<-mutate(hf, functional_seq = l75.jd -fopn.jd)
hf<-mutate(hf, inter_seq=bb.jd-fopn.jd)

##fopn and l75
by_species <- hf %>% group_by(species,year)
by_species<-by_species %>% summarise(mean(functional_seq, na.rm = TRUE), sd(functional_seq, na.rm = TRUE))

###fopn and lbb
by_species2 <- hf %>% group_by(species,year)
by_species2<-by_species2 %>% summarise(mean(inter_seq, na.rm = TRUE),sd(inter_seq, na.rm = TRUE))

##lbb and fbb
by_species3 <- hf %>% group_by(species,year)
by_species3<-by_species3 %>% summarise(mean(buds_seq, na.rm = TRUE),sd(buds_seq, na.rm = TRUE))

spec_tot<- hf %>% group_by(species)
spec_tot<-spec_tot %>% summarise(mean(functional_seq, na.rm = TRUE), sd(functional_seq, na.rm = TRUE))

spec_tot2<- hf %>% group_by(species)
spec_tot2<-spec_tot2 %>% summarise(mean(inter_seq, na.rm = TRUE), sd(functional_seq, na.rm = TRUE))

spec_tot3<- hf %>% group_by(species)
spec_tot3<-spec_tot3 %>% summarise(mean(buds_seq, na.rm = TRUE), sd(functional_seq, na.rm = TRUE))


#correlation between offset and overall growing degree days. ie can offselt of phenophases be explained by different growing degree needs?
#look at slopes



