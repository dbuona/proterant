#### Using HF species Where do flowering and leafing fall across species?

rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)
library(arm)
library(plyr)
library(broom)
library(tidyr)
setwd("~/Documents/git/proterant")
d<-read.csv("data/hf003-05-mean-ind.csv",header=TRUE)
d.spp<-read.csv("data/hf003-06-mean-spp.csv",header = TRUE)
###compture the overall mean flowering and leafing for each species
d<-filter(d, species %in% c( "ACPE","ACRU", "ACSA","BEAL","ILVE","NEMU","POTR","FRAM","QURU"))


sum<-group_by(d,species)
sum1<-dplyr::summarise(sum,ave.bb=mean(as.numeric(bb.jd),na.rm=TRUE))
sum2<-dplyr::summarise(sum,ave.fb=mean((fbb.jd),na.rm=TRUE))
?summarise
sum3<-dplyr::summarise(sum,ave.l75=mean(as.numeric(l75.jd),na.rm=TRUE))
sum4<-dplyr::summarise(sum,ave.fopn=mean(as.numeric(fopn.jd)na.rm=TRUE))

bb<-full_join(sum1,sum2,by="species")
bb<-gather(bb,phenophase,DOY,2:3)
bb$offset<-NA
bb<-offset

ggplot(bb,aes(x=DOY, y=1))+geom_point(aes(shape=phenophase, color=species))+theme_bw()



q<-ggplot(bb,aes(x=DOY, y=1))+geom_point(aes(shape=phenophase))
q+facet_wrap("species")

