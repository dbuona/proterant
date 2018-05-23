####Use harvard forest to evaluate annual variability hysteranthy
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)
library("randomForest")

####c### This is the final analysis file for hysteranthy anaylsis on MTSV as of 3/28/18.
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)
library("randomForest")

d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
unique(d$species)

a<-filter(d, species %in% c("ACRU","ACSA","AMSP" ,"BEAL" ,"BELE" ,"BEPA", "BEPO","FAGR","FRAM","POTR","PRSE","QUAL" ,"QURU" ,"QUVE","NYSY"))

#unique(hys$species) ##15 hysteranthous species (pro and syn)

a<-gather(a, phase,DOY,4:7)
a<-filter(a, phase %in% c("fopn.jd","l75.jd"))
ggplot(a,aes(year,DOY))+geom_point(aes(color=phase))+facet_wrap(~species)

hys<-filter(d
hys$offset<-NA
hys$offset<-hys$l75.jd-hys$fopn.jd
hys<-filter(hys, year<=2001)

ACRU<-filter(hys, species=="ACRU")
ACSA<-filter(hys, species=="ACSA")
BEAL<-filter(hys, species=="BEAL")
FRAM<-filter(hys, species=="FRAM")
QURU<-filter(hys, species=="QURU")
QUVE<-filter(hys, species=="QUVE")
POTR<-filter(hys, species=="POTR")
BEPO<-filter(hys, species=="BEPO")
QUER<-filter(hys, species %in% c("QUVE","QURU"))
ggplot(ACRU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
#ggplot(ACSA,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
#ggplot(BEAL,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
#ggplot(FRAM,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(QURU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(POTR,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(BEPO,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(QUVE,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
### Best
var<-filter(hys, species %in% c("QURU","ACRU","BEPO","POTR"))
var<-filter(var,!is.na(offset))
var<- var %>% group_by(tree.id) %>% summarise(sd(offset))

temp<-select(hys,tree.id,species)
temp<-unique(temp)
var.dat.<-left_join(var,temp)
colnames(var.dat.)<-c("tree.id","sd","species")

ggplot(var.dat., aes(species, sd))+stat_summary()
#year
ggplot(QUER,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd), color="pink")+geom_point(aes(year,bb.jd), color="green")+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd))+facet_wrap(~tree.id)
ggplot(ACRU,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd), color="pink")+geom_point(aes(year,bb.jd), color="green")+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd))+facet_wrap(~tree.id)
ggplot(BEPO,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd), color="pink")+geom_point(aes(year,bb.jd), color="green")+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd))+facet_wrap(~tree.id)
ggplot(QUVU,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd), color="pink")+geom_point(aes(year,bb.jd), color="green")+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd))+facet_wrap(~tree.id)
