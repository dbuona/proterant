rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/")

library(ape)
library(phytools)
library(brms)
library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")
library("jpeg")

d<-read.csv("Data/rosaceae.csv")

d<-filter(d,genus %in% c("Prunus","Rhododendron"))
colnames(d)


d$mean_size<-(d$flower_size_high+d$flos_per_low)/2
d$mean_num<-(d$flos_per_high+d$flos_per_low)/2
d$display<-d$mean_num*d$mean_size

a<-ggplot(d,aes(as.factor(hysteranthy),mean_size))+stat_summary()+
ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")+
  scale_x_discrete(labels=c("seranthous", "hysteranthous"))+facet_wrap(~genus,scale="free_y")  

b<-ggplot(d,aes(as.factor(hysteranthy),mean_num))+stat_summary()+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flowers/inflorescence")+
scale_x_discrete(labels=c("seranthous", "hysteranthous"))+facet_wrap(~genus,scale="free_y")  


c<-ggplot(d,aes(as.factor(hys.2),mean_size))+stat_summary()+ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")+facet_wrap(~genus) 
 

cc<-ggplot(d,aes(as.factor(hys.2),mean_num))+stat_summary()+ggthemes::theme_base()+xlab("FLS")+ylab("mean flowers/inflorescence")+
  scale_x_discrete(labels=c("seranthous", "synanthous", "hysteranthous", "hyst/syn"))  


table(d$hysteranthy)
table(d$hys.2)

d$hys.3<-d$hys.2
d$hys.3[which(d$hys.2=="1_0")]<-1
d$hys.3[which(d$hys.2=="0_-1")]<-1
d$hys.3[which(d$hys.2=="1_-1")]<-1
e<-ggplot(d,aes(as.factor(hys.3),mean_size))+stat_summary()+ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")+
  scale_x_discrete(labels=c("seranthous", "synanthous", "hysteranthous"))+facet_wrap(~genus) 

f<-ggplot(d,aes(as.factor(hys.3),mean_num))+stat_summary()+ggthemes::theme_base()+xlab("FLS")+ylab("mean flowers/inflorescence")+
  scale_x_discrete(labels=c("seranthous", "synanthous", "hysteranthous")) 

g<-ggplot(d,aes(as.factor(hysteranthy),display))+stat_summary()+
  ggthemes::theme_base()+xlab("FLS")+ylab("display volume (fl. size*fl. number")+
  scale_x_discrete(labels=c("seranthous", "hysteranthous"))+facet_wrap(~genus,scale="free_y")  


png(filename = "investment/prunusplot",width = 11, height=9,units = "in",res = 300)
ggpubr::ggarrange(a,b,g, nrow=3)
dev.off()
ggplot(d,aes(as.factor(hysteranthy),inflor_size_high))+stat_summary()+ggthemes::theme_base()+xlab("FLS")+ylab("inflorescence mid rib (mm)")+
  scale_x_discrete(labels=c("seranthous", "hysteranthous")) 


ggplot(d,aes(mean_num,mean_size))+geom_point()+facet_grid(genus~as.factor(hysteranthy))
