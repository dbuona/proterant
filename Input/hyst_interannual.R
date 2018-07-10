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



ACRU<-filter(d, species=="ACRU")
ACSA<-filter(d, species=="ACSA")
BEAL<-filter(d, species=="BEAL")
FRAM<-filter(d, species=="FRAM")
QURU<-filter(d, species=="QURU")
QUVE<-filter(d, species=="QUVE")
POTR<-filter(d, species=="POTR")
BEPO<-filter(d, species=="BEPO")
QUER<-filter(d, species %in% c("QUVE","QURU"))
ggplot(ACRU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(ACSA,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(BEAL,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(FRAM,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(QURU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(POTR,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(BEPO,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
ggplot(QUVE,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,l75.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+facet_wrap(~tree.id)
## do these for budf burst
ggplot(ACRU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red")+geom_point(aes(year,bb.jd), color="dark green")+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=bb.jd))+facet_wrap(~tree.id)
library("ggthemes")

QURU4<-filter(QURU,tree.id=="QURU-04")
ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fbb.jd), color="red",shape=8)+geom_point(aes(year,bb.jd), color="dark green",shape=18)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd))+theme_base()+ggtitle("Flowers bud burst with leaf budburst")+facet_wrap(~tree.id)
ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8)+geom_point(aes(year,bb.jd), color="dark green",shape=18)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=bb.jd))+theme_base()+ggtitle("Flowers open with leaf budburst")+facet_wrap(~tree.id)
ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8)+geom_point(aes(year,l75.jd), color="dark green",shape=18)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+theme_base()+ggtitle("Flowers open with leaf budburst")+facet_wrap(~tree.id)
??ggtitle
?shape()

##one tree only
##Fram-04
fram4<-filter(FRAM,tree.id=="FRAM-04")
ggplot(fram4,aes(year,fopn.jd))+geom_point(aes(year,fbb.jd), color="red",shape=8)+geom_point(aes(year,bb.jd), color="dark green",shape=18)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd))+theme_base()+ggtitle("Flowers bud burst with leaf budburst")+facet_wrap(~tree.id)
ggplot(fram4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8)+geom_point(aes(year,bb.jd), color="dark green",shape=18)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=bb.jd))+theme_base()+ggtitle("Flowers open with leaf budburst")+facet_wrap(~tree.id)
ggplot(fram4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8)+geom_point(aes(year,l75.jd), color="dark green",shape=18)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd))+theme_base()+ggtitle("Flowers open with leaf budburst")+facet_wrap(~tree.id)


ß### Best
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

###Individual interannual variability
d$offset<-d$fbb.jd-d$bb.jd
indies<-d %>% group_by(tree.id) %>% summarise(variation=sd(offset, na.rm=TRUE))
###withing year between individuals
years<-d%>% group_by(species,year) %>% summarise(physiological_variation = sd(offset,na.rm=TRUE)) 
d$offset2<-d$fopn.jd-d$l75.jd
indies2<-d %>% group_by(tree.id) %>% summarise(variation=sd(offset2, na.rm=TRUE))
years2<-d%>% group_by(species,year) %>% summarise(functional_variation = sd(offset2,na.rm=TRUE)) 

pop<-cbind(years,years2)
