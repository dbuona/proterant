

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
d$name[d$species=="ACPE"]<-"A. pensylvanicum"
d$name[d$species=="ACRU"]<-"A. rubrum"
d$name[d$species=="BEAL"]<-"B. alleghaniensis"
d$name[d$species=="FRAM"]<-"F. americana"
d$name[d$species=="NYSY"]<-"N. sylvatica"
d$name[d$species=="QURU"]<-"Q. rubra"


a<-filter(d, species %in% c("ACRU","BEAL","FRAM" ,"QURU","ACPE","NYSY" ))
colnames(a)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion (75%)","flower budburst","flower open","name")
unique(hys$species) ##15 hysteranthous species (pro and syn)

a<-gather(a, phase,DOY,4:7)

a$name <- factor(a$name, levels = c("A. rubrum","B. alleghaniensis","F. americana" ,"Q. rubra","A. pensylvanicum","N. sylvatica"))
pd<-position_dodge(width=0.0)

jpeg("..//figure/HFmeans.jpeg")
ggplot(a,(aes(name,DOY)))+stat_summary(fun.data = "mean_cl_boot",aes(color=phase,shape=phase),position=pd,size=0.78)+scale_color_manual(values=c("deeppink","firebrick1","lawngreen","darkgreen"))+scale_shape_manual(values=c(0,8,16,23))+theme_bw()+ylab("Day of Year")+xlab(NULL)#+theme(axis.text.x = element_text(angle = 300,hjust=0.5))
dev.off()



#a<-filter(a, phase %in% c("fopn.jd","l75.jd"))
#bet<-filter(a, species %in% c("BEAL"))
#ggplot(bet,aes(year,DOY))+geom_point(aes(shape=phase))+geom_smooth(method='lm',aes(,color=phase))+theme_base()
#ggplot(bet,aes(year,DOY))+geom_point(aes(shape=phase))+geom_smooth(method='lm',aes(,color=phase))+facet_wrap(~species)

###Alll Quercus and Fraxinus indivudals for proposal
QURU<-filter(d,species %in%c("QURU"))

QURU$budding<-QURU$fbb.jd-QURU$bb.jd
QURU<-filter(QURU,tree.id!="QURU-02")
QURU$FLS<-ifelse(QURU$budding<0,"hysteranthous","seranthous")
#QURU<-drop_na(QURU)
pd<-position_dodge(0.6)
a<-ggplot(QURU,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd,color=tree.id,group = row.names(QURU)),shape=8, size=4,position=pd)+geom_point(aes(year,bb.jd,color=tree.id,group = row.names(QURU)),shape=18,size=3,position=pd)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd)+theme_bw()+labs(y = "Day of year",color= "Tree I.D.")

QURU$funing<-QURU$fopn.jd-QURU$l75.jd

QURU$FLS<-ifelse(QURU$funing<0,"hysteranthous","seranthous")
#QURU<-drop_na(QURU)
pd<-position_dodge(0.6)
b<-ggplot(QURU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd,color=tree.id,group = row.names(QURU)),shape=8, size=4,position=pd)+geom_point(aes(year,l75.jd,color=tree.id,group = row.names(QURU)),shape=18,size=3,position=pd)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=l75.jd, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd)+theme_bw()+labs(y = "Day of year")

QURU$floping<-QURU$fopn.jd-QURU$bb.jd

QURU$FLS<-ifelse(QURU$floping<0,"hysteranthous","seranthous")
#QURU<-drop_na(QURU)
pd<-position_dodge(0.6)
c<-ggplot(QURU,aes(year,fopn.jd))+geom_point(aes(year,fopn.jd,color=tree.id,group = row.names(QURU)),shape=8, size=4,position=pd)+geom_point(aes(year,bb.jd,color=tree.id,group = row.names(QURU)),shape=18,size=3,position=pd)+geom_linerange(aes(x=year,ymin=fopn.jd,ymax=bb.jd, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd)+theme_bw()+labs(y = "Day of year")

grid.arrange(a,b,c)
