##basic exploration of harvard forest data
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
hf<-mutate(hf, buds_seq = bb.jd -fbb.jd)
hf<-mutate(hf, functional_seq = l75.jd -fopn.jd)
hf<-mutate(hf, inter_seq=bb.jd-fopn.jd)
library(dplyr)
by_species <- hf %>% group_by(species)
by_species<-by_species %>% summarise(mean(functional_seq, na.rm = TRUE))

by_species2 <- hf %>% group_by(species)
by_species2<-by_species2 %>% summarise(mean(inter_seq, na.rm = TRUE))

by_species3 <- hf %>% group_by(species)
by_species3<-by_species3 %>% summarise(mean(buds_seq, na.rm = TRUE))
