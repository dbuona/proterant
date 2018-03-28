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
#d<-filter(d, species %in% c( "ACPE","ACRU", "ACSA","BEAL","ILVE","NEMU","POTR","FRAM","QURU"))


sum<-group_by(d.spp,species)
#sum1<-dplyr::summarise(sum,ave.bb=mean(as.numeric(bb.jd),na.rm=TRUE))
#sum2<-dplyr::summarise(sum,ave.fb=mean((fbb.jd),na.rm=TRUE))

#sum3<-dplyr::summarise(sum,ave.l75=mean(as.numeric(l75.jd),na.rm=TRUE))
#sum4<-dplyr::summarise(sum,ave.fopn=mean(as.numeric(fopn.jd)na.rm=TRUE))

#bb<-full_join(sum1,sum2,by="species")
#bb<-gather(bb,phenophase,DOY,2:3)
bb<-gather(sum,phenophase,DOY,3:6)
bb<-filter(bb,phenophase==c("l75.jd","fopn.jd"))
ggplot(bb,aes(species,DOY))+stat_summary(aes(color=phenophase))

unique(bb$species)
trees<-c("ACPE","ACRU","ACSA","BEAL","BELE","BEPA","BEPO","CADE","FAGR","FRAM","NYSY","POTR","PRSE","QUAL","QURU","QUVE")
d.trees<-filter(bb, species %in% trees)
plot2<-ggplot(d.trees,aes(species,DOY))+stat_summary(aes(color=phenophase))+geom_abline(slope=0,intercept=142,color="red")+geom_abline(slope=0,intercept=148,color="blue")+ggtitle("Trees:Leaf and Flower Open")

lav<-as.data.frame(d.spp$l75.jd)
lav<-filter(lav,!is.na(d.spp$l75.jd))
quantile(lav$`d.spp$l75.jd`)
shrubs<-c("AMSP", "ARSP" , "COAL", "CRSP" , "HAVI" ,"ILVE" ,"KAAN", "KALA" ,"LYLI" ,"NEMU" ,"RHSP", "SAPU" , "VACO", "VIAL" ,"VICA")
d.shrub<-filter(bb, species %in% shrubs)
ggplot(d.shrub,aes(species,DOY))+stat_summary(aes(color=phenophase))

#####bud burst
sum<-group_by(d.spp,species)
names(sum)[names(sum) == 'bb.jd'] <- 'lbb.jd'
bb<-gather(sum,phenophase,DOY,3:6)
bb2<-filter(bb,phenophase==c("lbb.jd","fbb.jd"))
d.trees2<-filter(bb2, species %in% trees)
plot1<-ggplot(d.trees2,aes(species,DOY))+stat_summary(aes(color=phenophase))+geom_abline(slope=0,intercept=119,color="red")+geom_abline(slope=0,intercept=125,color="blue")+ggtitle("Trees:bb and FBB")

bav<-as.data.frame(d.spp$bb.jd)
bav<-filter(bav,!is.na(d.spp$bb.jd))
quantile(bav$`d.spp$bb.jd`)
median(bav$`d.spp$bb.jd`)

library(gridExtra)
grid.arrange(plot1,plot2, ncol=2)
