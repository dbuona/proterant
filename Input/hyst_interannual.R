

####c### This is the final analysis file for hysteranthy anaylsis on MTSV as of 3/28/18.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")

library("tidyverse")
library("ggplot2")
library(ggthemes)
library(gridExtra)

d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
unique(d$species)

a<-filter(d, species %in% c("ACRU","ACSA","AMSP" ,"BEAL" ,"BELE" ,"BEPA", "BEPO","FRAM","POTR" ,"QURU" ,"QUVE"))

#unique(hys$species) ##15 hysteranthous species (pro and syn)

#a<-gather(a, phase,DOY,4:7)
#a<-filter(a, phase %in% c("fopn.jd","l75.jd"))
#bet<-filter(a, species %in% c("BEAL"))
#ggplot(bet,aes(year,DOY))+geom_point(aes(shape=phase))+geom_smooth(method='lm',aes(,color=phase))+theme_base()
#ggplot(bet,aes(year,DOY))+geom_point(aes(shape=phase))+geom_smooth(method='lm',aes(,color=phase))+facet_wrap(~species)


QURU4<-filter(d,tree.id %in%c("FRAM-04"))
QURU4$budding<-QURU4$fbb.jd-QURU4$bb.jd
QURU4$bud2<-ifelse(QURU4$budding<0,"hysteranthous","seranthous")
b<-ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fbb.jd), color="red",shape=10, size=1.5)+geom_point(aes(year,bb.jd), color="dark green",shape=1,size=3)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd, linetype=bud2))+theme_base()+facet_wrap(~tree.id)+theme(legend.position="none")

QURU4$flobud<-QURU4$fopn.jd-QURU4$bb.jd
QURU4$bud3<-ifelse(QURU4$flobud<0,"hysteranthous","seranthous")
c<-ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8, size=2)+geom_point(aes(year,bb.jd), color="dark green",shape=1,size=3)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=bb.jd, linetype=bud3))+theme_base()+facet_wrap(~tree.id)+theme(legend.position="none")

QURU4$dev<-QURU4$fopn.jd-QURU4$l75.jd
QURU4$bud5<-ifelse(QURU4$dev<0,"hysteranthous","seranthous")
e<-ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8, size=2)+geom_point(aes(year,l75.jd), color="dark green",shape=5,size=3)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd, linetype=bud5))+theme_base()+facet_wrap(~tree.id)+theme(legend.position="none")

y<-grid.arrange(b,c,e,ncol=3)

QURU4<-filter(d,tree.id %in%c("QURU-04"))
QURU4$budding<-QURU4$fbb.jd-QURU4$bb.jd
QURU4$bud2<-ifelse(QURU4$budding<0,"hysteranthous","seranthous")
f<-ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fbb.jd), color="red",shape=10, size=1.5)+geom_point(aes(year,bb.jd), color="dark green",shape=1,size=3)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd, linetype=bud2))+theme_base()+facet_wrap(~tree.id)+theme(legend.position="none")

QURU4$flobud<-QURU4$fopn.jd-QURU4$bb.jd
QURU4$bud3<-ifelse(QURU4$flobud<0,"hysteranthous","seranthous")
g<-ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8, size=2)+geom_point(aes(year,bb.jd), color="dark green",shape=1,size=3)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=bb.jd, linetype=bud3))+theme_base()+facet_wrap(~tree.id)+theme(legend.position="none")

QURU4$dev<-QURU4$fopn.jd-QURU4$l75.jd
QURU4$bud5<-ifelse(QURU4$dev<0,"hysteranthous","seranthous")
h<-ggplot(QURU4,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd), color="red",shape=8, size=2)+geom_point(aes(year,l75.jd), color="dark green",shape=5,size=3)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd, linetype=bud5))+theme_base()+facet_wrap(~tree.id)+theme(legend.position="none")
x<-grid.arrange(f,g,h,ncol=3)



