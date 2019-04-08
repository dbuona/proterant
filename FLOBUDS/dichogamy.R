###this model explores changes in dichogamy in a subset of species

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(tidyverse)
library(brms)
library(rstan)
library(arm)
library(rstanarm)
library(tibble)
library(ggstance)

setwd("~/Documents/git/proterant/FLOBUDS")

d<-read.csv("dicogamy.csv",header=TRUE)

d<-unite(d,treatment,Force,Light,Chill,sep="",remove=FALSE)
d$Light<-ifelse(d$Light=="L","xL","S")
###name treatments numeric/continuous
d$photoperiod<-ifelse(d$Light=="xL",12,8)
d$temp_day<-ifelse(d$Force=="W",24,18)
d$temp_night<-ifelse(d$Force=="W",18,12)
d$chilldays<-ifelse(d$Chill==0,28,56)

###center predictors
d$p_cent<-d$photoperiod/mean(d$photoperiod)
d$f_cent<-d$temp_day/mean(d$temp_day)
d$c_cent<-d$chilldays/mean(d$chilldays)

d$p.z<-d$photoperiod-mean(d$photoperiod)/sd(d$photoperiod)
d$f.z<-d$temp_day-mean(d$temp_day)/sd(d$temp_day)
d$c.z<-d$chilldays-mean(d$chilldays)/sd(d$chilldays)

goo<- d %>% group_by(GEN.SPA,treatment) %>% summarise(mean.gyn=mean(flo_dayF,na.rm=TRUE))
goo2<-d %>% group_by(GEN.SPA,treatment) %>% summarise(mean.anth=mean(flo_dayM,na.rm=TRUE))
googoo<-left_join(goo,goo2)
googoo$dichog<-googoo$mean.gyn-googoo$mean.anth

cor<-filter(d,GEN.SPA=="COR.COR")
com<-filter(d,GEN.SPA=="COM.PER")


ladylikelinood<-brm(gynlike~p.z+f.z+c.z,data=com, family=bernoulli(link="logit"))
summary(ladylikelinood)
brolikelinood<-brm(anthlike~p.z+f.z+c.z,data=com, family=bernoulli(link="logit"))


ladylikelinood.cor<-brm(gynlike~p.z+f.z+c.z+p.z:f.z+p.z:c.z+f.z:c.z,data=cor, family=bernoulli(link="logit"))
summary(ladylikelinood.cor)
brolikelinood.cor<-brm(anthlike~p.z+f.z+c.z,data=cor, family=bernoulli(link="logit"))
summary(brolikelinood.cor)


lady<-brm(flo_dayF~p.z+f.z+c.z,data=com)
bro<-brm(flo_dayM~~p.z+f.z,data=com) ### doesn't flower without chilling
summary(lady)
summary(bro)
lady2<-brm(flo_dayF~p.z+f.z+c.z,data=cor)
bro2<-brm(flo_dayM~~p.z+f.z+c.z,data=cor)


summary(lady2)
summary(bro2)
gynsums<- d %>% group_by(GEN.SPA,treatment) %>%  summarise(mean.gyn=mean(flo_dayF, na.rm=TRUE)) 
anthsums<- d %>% group_by(GEN.SPA,treatment)%>%  summarise(mean.anth=mean(flo_dayM,na.rm=TRUE))
gynsums2<- d %>% group_by(GEN.SPA,treatment) %>%  summarise(sd.gyn=sd(flo_dayF, na.rm=TRUE)) 
anthsums2<- d %>% group_by(GEN.SPA,treatment)%>%  summarise(sd.anth=sd(flo_dayM,na.rm=TRUE))
dicho<-left_join(gynsums,gynsums2)
dich2<-left_join(anthsums,anthsums2)
dicho<-left_join(dicho,dich2)
dicho$dicogamy<-dicho$mean.anth-dicho$mean.gyn

ggplot(dicho, aes(treatment, dicogamy))+geom_boxplot(aes(color=GEN.SPA))
dicho$dichogamybin<-ifelse(dicho$dicogamy<=0,1,0)
ggplot(dicho, aes(treatment,dichogamybin))+geom_point()+facet_wrap(~GEN.SPA)

summary(lady)

         