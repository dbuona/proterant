###TO DO run models without ACE.SAC and BET.ALL to see if figures work better
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library("brms")
library(rstan)
library(arm)
library(ggstance)
library(dplyr)
library(ggthemes)
library("Hmisc")
library(brms)
library(broom)
library(RColorBrewer)
library("tidybayes")

#options(device = "quartz")
setwd("~/Documents/git/proterant/FLOBUDS")
load("writing/flobud.main.mods.Rda")
dat<-read.csv("flobudsdata.use.csv",header = TRUE)

dat$Light<-ifelse(dat$Light=="S",0,1)
dat$Force<-ifelse(dat$Force=="C",0,1)

dat<-filter(dat, !GEN.SPA %in% c("AME.SPP","BET.SPP", "BET.ALL","ACE.SAC")) ## remove species with no flowering



#####plasticity models for 3 bbch stages 
bb.int<-get_prior(budburst.9.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())

mod.bb.int<-brm(budburst.9. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                   data = dat, family = gaussian(),
                   iter= 4000,
                   warmup = 3000)   



flo.int<-get_prior(flo_day~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())
mod.flo.int<-brm(flo_day~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3000)   




lo.int<-get_prior(leaf_day.15.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())

mod.lo.int<-brm(leaf_day.15. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3000)   


save.image("flobud.main.mods.Rda")

summary(mod.flo.int)
summary(mod.bb.int)
summary(mod.lo.int)
fixef(mod.lo.int)
fixef(mod.flo.int)

dat %>% dplyr::group_by(GEN.SPA) %>% dplyr::summarise(meandvr=mean(leaf_day.15.-budburst.9.,na.rm=TRUE))

sd(dat$leaf_day.15.-dat$budburst.9.,na.rm=TRUE)

bayestestR::effective_sample(mod.flo.int,effects = c("all"))
bayestestR::effective_sample(mod.bb.int)





####################################



###now plots for comparing file observation to mine
