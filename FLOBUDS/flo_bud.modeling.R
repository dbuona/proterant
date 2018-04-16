
rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(ggplot2)
library(brms)
library(rstan)
library(arm)
library(rstanarm)
library(tibble)

setwd("~/Documents/git/proterant/FLOBUDS")

d<-read.csv("first.flobud.data.csv",header=TRUE)
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

ggplot(d,aes(flo_day))+geom_density()

mod<-brm(cbind(leaf_day,flo_day) ~ Light+Chill+Force+Light:Chill+Light:Force+Force:Chill+(Light +Chill+Force+Light:Chill+Light:Force+Force:Chill|p|GEN.SPA),
               data = d, family = gaussian, 
               iter= 3000,
               warmup = 2000,
               cores = 4)

summary(mod)
coef(mod [1])



modA<-brm(flo_day ~ photoperiod +Chill+temp_day+photoperiod:Chill+photoperiod:temp_day+temp_day:Chill+( photoperiod +Chill+temp_day+photoperiod:Chill+photoperiod:temp_day+temp_day:Chill|GEN.SPA),
                   data = d, family = gaussian, 
                  iter= 3000,
                  warmup = 2000,
                  cores = 4)


modA<-brm(cbind(leaf_day,flo_day) ~ p_cent +c_cent+f_cent+p_cent:c_cent+p_cent:f_cent+f_cent:c_cent+(1+ p_cent +c_cent+f_cent+p_cent:c_cent+p_cent:f_cent+f_cent:c_cent|p|GEN.SPA),
          data = d, family = gaussian, 
          iter= 3000,
          warmup = 2500,
          cores = 4)
brms::pp_check(modA, resp= "leafday")
brms::pp_check(modA, resp= "floday") 
summary(modA)
ranef(modA)
coef(modA)
bayes_R2(modA)


