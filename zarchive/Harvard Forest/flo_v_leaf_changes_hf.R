###plotting differences in trends.
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
#hf<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
mean_hf<-read.csv("hf003-06-mean-spp.csv",header = TRUE)
dat<-gather(mean_hf,phenophase, eventday, bb.jd:fopn.jd)
dat<-dat %>% group_by(year, species,phenophase) %>% summarise(meaneventday=mean(eventday))
b <- ggplot(dat, aes(x = year, y = meaneventday, color=phenophase))+ geom_point()
b+geom_smooth(method = "lm")

###now indiv species
mean_hf<-filter(mean_hf, species %in% c( "ACPE","ACRU", "ACSA","BEAL","FRAM","QURU"))
dat<-gather(mean_hf,phenophase, eventday, bb.jd:fopn.jd)
dat<-dat %>% group_by(year, species,phenophase) %>% summarise(meaneventday=mean(eventday))
b <- ggplot(dat, aes(x = year, y = meaneventday, color=phenophase))+ geom_point()
b+geom_smooth(method = "lm")+facet_wrap(~species)
