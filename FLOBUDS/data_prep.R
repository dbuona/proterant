####data leadin 6 Nov 2018
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(tidyverse)
library("brms")
library(rstan)
library(arm)
library(rstanarm)
library(tibble)
library(ggstance)
library(survival)
library(sur)
library(survminer)
library(ggthemes)
library("Hmisc")

setwd("~/Documents/git/proterant/FLOBUDS")
#load("flobud_mods.RData")
d<-read.csv("input/first.event.dat.csv",header=TRUE)
###get rid of useless colummns
colnames(d)
d<-dplyr::select(d, -X.1)
d<-dplyr::select(d, -X)

#### make light differsent so variables  appear in order
#d$Light<-ifelse(d$Light=="L","XL","S")
###name treatments numeric/continuous
d$photoperiod<-ifelse(d$Light=="L",12,8)
d$temp_day<-ifelse(d$Force=="W",24,18)
d$temp_night<-ifelse(d$Force=="W",18,12)
d$chilldays<-ifelse(d$Chill==0,28,56)

###center predictors
d$p_scale<-d$photoperiod-mean(d$photoperiod)/2*(sd(d$photoperiod))
d$f_scale<-d$temp_day-mean(d$temp_day)/2*(sd(d$temp_day))
d$c_scale<-d$chilldays-mean(d$chilldays)/2*(sd(d$chilldays))

###what is the distrbution of response variables
ggplot(d,aes(flo_day))+geom_density()
ggplot(d,aes(leaf_day))+geom_density()
ggplot(d,aes(Lbb_day))+geom_density()
ggplot(d,aes(Lexpand_day))+geom_density()
d<-unite(d, treatment, Force,Light, Chill, sep= "",remove = FALSE)

########Basic plots for each phenophase#############################################

phased<-gather(d,phase,DOY,8:12)

bbview<-filter(phased, phase %in% c("Lbb_day","flo_day"))
expview<-filter(phased, phase %in% c("Lexpand_day","flo_day"))
leafview<-filter(phased, phase %in% c("leaf_day","flo_day"))

###leafout and flowering
bigsp<-filter(leafview, !GEN.SPA %in% c("AME.SPP","BET.SPP"))
unique(bigsp$GEN.SPA)

library(ggstance)
pd=position_dodge(0.1)
bigsp$phase[bigsp$phase=="flo_day"]<-"flower"
bigsp$phase[bigsp$phase=="leaf_day"]<-"leaf"             
bigsp$Light[bigsp$Light=="XL"]<-"long photoperiod"
bigsp$Light[bigsp$Light=="S"]<-"short photoperiod"
bigsp$Chill[bigsp$Chill=="0"]<-"short chilling"
bigsp$Chill[bigsp$Chill=="1"]<-"long chilling"
bigsp$Force[bigsp$Force=="C"]<-"low forcing"
bigsp$Force[bigsp$Force=="W"]<-"high forcing"

p<-ggplot(bigsp,aes(GEN.SPA,as.numeric(DOY)))+geom_point(aes(shape=phase,color=phase),size=0.8)+ylab("days to event")+xlab("species")+stat_summary(fun.data = "mean_cl_boot",aes(shape=phase,color=phase),position=pd)+facet_grid(Force~Light~Chill)+theme_bw()
pp<-p+theme(axis.text.x = element_text(size=8,angle = 300, hjust = 0))

###################survival analysis###########Kaplan-Meier########################
viv<-filter(d,Dead.alive %in% c("A","?"))
table(viv$treatment)

table(viv$treatment)
table(d$treatment) 

vivo.full<-gather(viv,phase,DOY,9:12)

###do it for bud burst
vivo<-filter(vivo.full, phase %in% c("Lbb_day","flo_day"))
vivo$DOY<-ifelse(is.na(vivo$DOY),120,vivo$DOY)
vivo$surv<-ifelse(vivo$DOY==120,1,0)

vivo$floposs<-ifelse(vivo$Flo.poss.=="N",0,1)

###This is only twigs where we deemed flowering to even bee possible
vivo2<-filter(vivo,floposs==1)

###vivo2 is main datasheet for buds and flower
###chayim 2 is for leaf out and flowering
unique(vivo.full$phase)
hayim<-dplyr::filter(vivo.full, phase %in% c("leaf_day","flo_day"))
hayim$DOY<-ifelse(is.na(hayim$DOY),120,hayim$DOY)
hayim$surv<-ifelse(hayim$DOY==120,1,0)
hayim$floposs<-ifelse(hayim$Flo.poss.=="N",0,1)
###This is only twigs where we deemed flowering to even bee possible
hayim2<-filter(hayim,floposs==1)
write.csv(vivo2,"budburst_survival_data")
write.csv(hayim2, "flo_exapand_survival_data")
