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
d<-read.csv("first.event.dat.csv",header=TRUE)
###get rid of useless colummns
colnames(d)
d<-dplyr::select(d, -X)

#### make light differsent so variables  appear in order
#d$Light<-ifelse(d$Light=="L","XL","S")
###name treatments numeric/continuous
d$photoperiod<-ifelse(d$Light=="L",12,8)
d$temp_day<-ifelse(d$Force=="W",24,18)
d$temp_night<-ifelse(d$Force=="W",18,12)
d$chilldays<-ifelse(d$Chill==0,28,56)

###zscore predictors
d$p_z<-d$photoperiod-mean(d$photoperiod)/(sd(d$photoperiod))
d$f_z<-d$temp_day-mean(d$temp_day)/(sd(d$temp_day))
d$c_z<-d$chilldays-mean(d$chilldays)/(sd(d$chilldays))

###what is the distrbution of response variables
ggplot(d,aes(flo_day.60.))+geom_density()
ggplot(d,aes(leaf_day.15.))+geom_density()
ggplot(d,aes(Lbb_day.9.))+geom_density()
ggplot(d,aes(Lexpand_day.11.))+geom_density()
d<-unite(d, treatment, Force,Light, Chill, sep= "",remove = FALSE)

########Basic plots for each phenophase#############################################

phased<-gather(d,phase,DOY,8:12)

bbview<-filter(phased, phase %in% c("Lbb_day.9.","flo_day.60."))
expview<-filter(phased, phase %in% c("Lexpand_day.11.","flo_day.60."))
leafview<-filter(phased, phase %in% c("leaf_day.15.","flo_day.60."))

###leafout and flowering
bigsp<-filter(leafview, !GEN.SPA %in% c("AME.SPP","BET.SPP"))
unique(bigsp$GEN.SPA)

library(ggstance)
pd=position_dodge(0.1)
bigsp$phase[bigsp$phase=="flo_day.60."]<-"flower"
bigsp$phase[bigsp$phase=="leaf_day.15."]<-"leaf"             
bigsp$Light[bigsp$Light=="XL"]<-"long photoperiod"
bigsp$Light[bigsp$Light=="S"]<-"short photoperiod"
bigsp$Chill[bigsp$Chill=="0"]<-"short chilling"
bigsp$Chill[bigsp$Chill=="1"]<-"long chilling"
bigsp$Force[bigsp$Force=="C"]<-"low forcing"
bigsp$Force[bigsp$Force=="W"]<-"high forcing"

pd=position_dodgev(height=0.4)
p<-ggplot(bigsp,aes(GEN.SPA,as.numeric(DOY)))+geom_point(aes(shape=phase,color=phase),size=0.8)+ylab("days to event")+xlab("species")+stat_summary(fun.data = "mean_cl_boot",aes(shape=phase,color=phase),position=pd)+facet_grid(Force~Light~Chill)+theme_bw()
pp<-p+theme(axis.text.x = element_text(size=8,angle = 300, hjust = 0))

###################survival analysis###########Kaplan-Meier########################
viv<-filter(d,Dead.alive %in% c("A","?"))
table(viv$treatment)

table(viv$treatment)
table(d$treatment) 

vivo.full<-gather(viv,phase,DOY,9:12)

###do it for bud burst
vivo<-filter(vivo.full, phase %in% c("Lbb_day.9.","flo_day.60."))
vivo$DOY<-ifelse(is.na(vivo$DOY),120,vivo$DOY)
vivo$surv<-ifelse(vivo$DOY==120,1,0)

vivo$floposs<-ifelse(vivo$Flo.poss.=="N",0,1)

###This is only twigs where we deemed flowering to even bee possible
vivo2<-filter(vivo,floposs==1)

###vivo2 is main datasheet for buds and flower
###chayim 2 is for leaf out and flowering

hayim<-dplyr::filter(vivo.full, phase %in% c("Lexpand_day.11.","flo_day.60."))
hayim$DOY<-ifelse(is.na(hayim$DOY),120,hayim$DOY)
hayim$surv<-ifelse(hayim$DOY==120,1,0)
hayim$floposs<-ifelse(hayim$Flo.poss.=="N",0,1)
###This is only twigs where we deemed flowering to even bee possible
hayim2<-filter(hayim,floposs==1)
unique(vivo.full$phase)
vita<-dplyr::filter(vivo.full, phase %in% c("leaf_day.9.","flo_day.60."))
vita$DOY<-ifelse(is.na(vita$DOY),120,vita$DOY)
vita$surv<-ifelse(vita$DOY==120,1,0)
vita$floposs<-ifelse(vita$Flo.poss.=="N",0,1)
###This is only twigs where we deemed flowering to even bee possible
vita2<-filter(vita,floposs==1)



write.csv(vivo2,"budburst_survival_data.csv")
write.csv(hayim2, "flo_exapand_survival_data.csv")
write.csv(vita2,"flo_leafout_survival_data.csv")
