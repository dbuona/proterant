###TO DO run models without ACE.SAC and BET.ALL to see if figures work better
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
library(brms)
library(broom)
library(RColorBrewer)

#options(device = "quartz")
setwd("~/Documents/git/proterant/FLOBUDS")
load("new_flobud.mods.Rda")
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


###summarize main effects of these models
library(xtable)
xtable(fixef(mod.bb.int,probs = c(.25,.75)))
xtable(fixef(mod.flo.int,probs = c(.25,.75)),caption = "goo",label ="goo" )

xtable(modoutput1)
xtable(modoutput2)
 ##sp effect



########## THis make a plot#################
figpath <- "Plots"

cols <- adjustcolor("indianred3", alpha.f = 0.3) 
my.pal <- rep(brewer.pal(n = 10, name = "Paired"), 8)
# display.brewer.all()
alphahere = 1

xlab <- "Model estimate of change in phenophase day"

spp <- unique(dat$GEN.SPA)

modelhere <-mod.bb.int
modelhere2<-mod.flo.int
leg.txt <- c("reproductive","vegetative")
source("exp_muplot_brms.R")
source("prep4plot.R")



muplotfx(modelhere,modelhere2, "budburstvsflowering", 8, 8, c(0,6), c(-50, 30) , 40, 3.5,40,4.5)
dev.off()
####################################



###now plots for comparing file observation to mine
small<-filter(dat,!is.na(flo_day))
small<-filter(small,!is.na(budburst.9.))

small$FLS<-small$budburst.9.-small$flo_day

small.fls<-get_prior(FLS~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = small, family = gaussian())
mod.FLS.small<-brm(FLS ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                  data = small, family = gaussian(),control = list(adapt_delta = 0.95),
                  iter= 7000,
                  warmup = 6000)



####climate change projections
### new data for predicting climat change
new.data<-data.frame(GEN.SPA=rep(unique(dat$GEN.SPA),9),
                     Force=rep(c(0,1),each=45),
                     Chill=rep(c(.67,1,0),30),
                     Light=rep(c(1),90))


#####now do flowering a leafing seperately
predictionflo<-predict(mod.flo.int,newdata=new.data,probs = c(.055,.945,.25,.75))
predy<-cbind(new.data,predictionflo)
predy$scenario<-NA
predy$scenario[which(predy$Force==0 & predy$Chill==.67)]<-"historic"
predy$scenario[which(predy$Force==1 & predy$Chill==.67)]<-"warm 5"
predy$scenario[which(predy$Force==1 & predy$Chill== 1)]<-"5+chill"
predy$scenario[which(predy$Force==1 & predy$Chill== 0)]<-"5-chill"

unique(predy$scenario)
predy <- na.omit(predy)
predy %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5", "5-chill", "5+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes())+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0))+facet_wrap(~GEN.SPA)+theme_bw()+geom_hline(yintercept=0)



predictionbb<-predict(mod.bb.int,newdata=new.data,probs = c(.055,.945,.25,.75))
predy2<-cbind(new.data,predictionbb)
predy2$scenario<-NA
predy2$scenario[which(predy2$Force==0 & predy2$Chill==.67)]<-"historic"
predy2$scenario[which(predy2$Force==1 & predy2$Chill==.67)]<-"warm 5"
predy2$scenario[which(predy2$Force==1 & predy2$Chill== 1)]<-"5+chill"
predy2$scenario[which(predy2$Force==1 & predy2$Chill== 0)]<-"5-chill"

unique(predy2$scenario)
predy2 <- na.omit(predy2)

predy$phase<-"flower"
predy2$phase<-"budburst"
predybig<-rbind(predy,predy2)



setEPS()
postscript("Plots/climpredictions.eps",width = 12, height = 8)
predybig %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=phase,shape=phase),size=2.5)+geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0,color=phase),linetype="solid")+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=phase),linetype="dotted")+facet_wrap(~GEN.SPA,scale="free",ncol=3)+
  ggthemes::theme_base(base_size = 10)+scale_color_manual(values=c("darkgreen","darkorchid3"))+scale_shape_manual(values=c(17,19))
dev.off()

jpeg("Plots/climpredictions_fixedscale.jpeg",width = 6, height = 10, units = 'in', res=300)
predybig %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=phase,shape=phase),size=2.5)+geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0,color=phase),linetype="solid",alpha=0.8)+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=phase),linetype="dotted",alpha=0.8)+facet_wrap(~GEN.SPA,scale="fixed",ncol=2)+theme_bw()+scale_color_manual(values=c("darkgreen","hotpink"))
dev.off()



###prediction plots
HFreal<-read.csv(file = "..//Data/hf003-05-mean-ind.csv")
HFreal$FLS<-HFreal$bb.jd-HFreal$fopn.jd
HFrealmeans<-HFreal %>% group_by(species) %>% dplyr::summarize(Estimate=mean(FLS,na.rm = TRUE),Q94.5=max(FLS,na.rm = TRUE),Q5.5=min(FLS,na.rm = TRUE))
HFrealmeans<-filter(HFrealmeans,species %in% c('ACPE',"ACRU","NEMU","ILVE","VACO"))
HFrealmeans$GEN.SPA<-NA

HFrealmeans$GEN.SPA[which(HFrealmeans$species=="ACPE")]<-"ACE.PEN"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="ACRU")]<-"ACE.RUB"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="NEMU")]<-"ILE.MUC"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="ILVE")]<-"ILE.VER"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="VACO")]<-"VAC.COR"
HFrealmeans$scenario<-"field"
HFrealmeans<-dplyr::select(HFrealmeans,-species)
coef(mod.FLS.small)
new.data<-data.frame(GEN.SPA=rep(unique(small$GEN.SPA),9),
                     Force=rep(c(0,1),each=45),
                     Chill=rep(c(.67,1,0),30),
                     Light=rep(c(1),90))



prediction<-predict(mod.FLS.small,newdata=new.data,probs = c(.055,.945))
predy<-cbind(new.data,prediction)


predy$scenario<-NA
predy$scenario[which(predy$Force==0 & predy$Chill==.67)]<-"historic"
predy$scenario[which(predy$Force==1 & predy$Chill==.67)]<-"warm 5"
predy$scenario[which(predy$Force==1 & predy$Chill== 1)]<-"5+chill"
predy$scenario[which(predy$Force==1 & predy$Chill== 0)]<-"5-chill"
unique(predy$scenario)
predy <- na.omit(predy) 
predy<-dplyr::select(predy,GEN.SPA,Estimate,Q5.5,Q94.5,scenario)

predy<-rbind(predy,HFrealmeans)

predy$grouper<-NA
predy$grouper[which(predy$scenario %in% c("field"))]<-"historic-observed"
predy$grouper[which(predy$scenario %in% c("historic"))]<-"historic-predicted"
predy$grouper[which(predy$scenario %in% c("warm 5", "warm 10"))]<-"warm only"
predy$grouper[which(predy$scenario %in% c("5-chill", "10-chill"))]<-"warm,reduce chill"
predy$grouper[which(predy$scenario %in% c("5+chill", "10+chill"))]<-"warm,increase chill"

predy<-filter(predy,GEN.SPA %in% c("ACE.PEN","ACE.RUB","ILE.MUC","ILE.VER","VAC.COR"))
predy %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("field","historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=grouper))+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=grouper))+facet_wrap(~GEN.SPA)+theme_bw()+geom_hline(yintercept=0)

prey<-filter(predy,scenario %in% c("field","historic"))
prey$scenario[which(prey$scenario %in% c("field"))]<-"observed"
prey$scenario[which(prey$scenario %in% c("historic"))]<-"predicted"

jpeg("Plots/fieldmodcomparisions_freescale.jpeg",width = 5, height = 6, units = 'in', res=300)
ggplot(prey,aes(scenario,Estimate))+geom_point(aes())+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0))+facet_wrap(~GEN.SPA,scales="free_y")+theme_bw()+geom_hline(yintercept=0)
dev.off()
jpeg("Plots/fieldmodcomparisions.jpeg",width = 5, height = 6, units = 'in', res=300)
ggplot(prey,aes(scenario,Estimate))+geom_point(aes())+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0))+facet_wrap(~GEN.SPA)+theme_bw()+geom_hline(yintercept=0)
dev.off()

save.image("new_flobud.mods.Rda")