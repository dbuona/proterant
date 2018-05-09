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
?unite()
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

count(d,!is.na(flo_dayF)) #31
count(d,!is.na(flo_dayM)) #36

cor<-filter(d,GEN.SPA=="COR.COR")
count(cor,!is.na(flo_dayF)) #12
count(cor,!is.na(flo_dayM)) #15

com<-filter(d,GEN.SPA=="COM.PER")
count(com,!is.na(flo_dayF)) #19
count(com,!is.na(flo_dayM)) #21

lady<-lm(flo_dayF~Light+Chill+Force+Light:Chill+Light:Force+Force:Chill,data=cor)

bro<-lm(flo_dayM~Light+Chill+Force+Light:Chill+Light:Force+Force:Chill,data=cor)

summary(lady)
summary(bro)

         