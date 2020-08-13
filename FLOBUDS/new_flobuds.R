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

dat<-filter(dat, !GEN.SPA %in% c("AME.SPP","BET.SPP", "BET.ALL","ACE.SAC"))

rawplot<-gather(dat,phase,DOY,3:4)
rawplot$Forcing<-ifelse(rawplot$Force==1,"high forcing","low forcing")
rawplot$Photoperiod<-ifelse(rawplot$Light==1,"long photoperiod","short photoperiod")
rawplot$Chilling<-ifelse(rawplot$Chill==1,"long chilling","short chilling")
rawplot$phenophase<-ifelse(rawplot$phase=="budburst.9.","budburst","flowering")

#setEPS()
#postscript("Plots/rawdataplots.eps",width = 10, height = 10)
jpeg("Plots/rawdataplots.jpg",width = 10, height = 10,units='in',res=200)
ggplot(rawplot,aes(GEN.SPA,DOY))+stat_summary(aes(color=phenophase,shape=phenophase))+geom_point(aes(color=phenophase,shape=phenophase),size=1,alpha=0.6)+facet_grid(Forcing~Photoperiod~Chilling,scales = "free_y")+
  scale_color_manual(values=c("darkgreen","darkorchid3"))+scale_shape_manual(values=c(17,19))+  
  theme_bw(base_size = 11,base_line_size = .2)+ylab("Day of Experiment")+theme(axis.text.x = element_text(angle = 30,hjust = 0.9))+xlab("")
dev.off()

vacplot<-filter(rawplot,GEN.SPA=="VAC.COR")
ggplot(vacplot,aes(Chilling,DOY))+stat_summary(aes(color=phenophase,shape=phenophase),size=.6)+geom_point(aes(color=phenophase,shape=phenophase),size=1,alpha=0.6)+facet_grid(Photoperiod~Forcing)+
  theme_bw(base_size = 10,base_line_size = .2)+scale_color_manual(values=c("darkgreen","darkorchid3"))+scale_shape_manual(values=c(17,19))+
  ylab("Day of Experiment")+xlab("V.corymbosum")+theme(axis.title.x = element_text(face = "italic"))

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


###check again only complete cases
small<-filter(dat,!is.na(flo_day))
small<-filter(small,!is.na(budburst.9.))
bb.small<-get_prior(budburst.9.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = small, family = gaussian())

mod.bb.small<-brm(budburst.9. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = small, family = gaussian(),control = list(adapt_delta = 0.95),
                iter= 7000,
                warmup = 6000)

fo.small<-get_prior(flo_day~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = small, family = gaussian())
mod.fo.small<-brm(flo_day ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                  data = small, family = gaussian(),control = list(adapt_delta = 0.95),
                  iter= 7000,
                  warmup = 6000)

summary(mod.bb.small)

lo.int<-get_prior(leaf_day.15.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())
mod.lo.int<-brm(leaf_day.15. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3000)   

dev.new()
pp_check(mod.bb.int)
##########
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

muplotfx(modelhere,modelhere2, "Intecept-budburst vs. flowering", 8, 8, c(6,8), c(30, 130) , 40, 3.5,40,4.5)
dev.off()

modelhere <-mod.lo.int
source("prep4plot.R")
source("exp_muplot_brms.R")

muplotfx(modelhere,modelhere2, "leafout vs. flowering", 8, 8, c(0,7), c(-50, 120) , 130, 3.5)
dev.off()

modelhere <-mod.bb.int
modelhere2<-mod.bb.small
source("prep4plot.R")
source("exp_muplot_brms.R")
muplotfx(modelhere,modelhere2, "budburst pool vs. complete cases", 8, 8, c(0,7), c(-50, 120) , 130, 3.5)
dev.off()

modelhere <-mod.flo.int
modelhere2<-mod.fo.small
source("prep4plot.R")
source("exp_muplot_brms.R")
muplotfx(modelhere,modelhere2, "flowering pool vs. complete cases", 8, 8, c(0,7), c(-50, 120) , 130, 3.5)
dev.off()

small$FLS<-small$budburst.9.-small$flo_day
small$FLS2<-small$leaf_day.15.-small$flo_day
small.fls<-get_prior(FLS~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = small, family = gaussian())
mod.FLS.small<-brm(FLS ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                  data = small, family = gaussian(),control = list(adapt_delta = 0.95),
                  iter= 7000,
                  warmup = 6000)

small.fls2<-get_prior(FLS2~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = small, family = gaussian())
mod.FLS2.small<-brm(FLS2 ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                   data = small, family = gaussian(),control = list(adapt_delta = 0.95),
                   iter= 7000,
                   warmup = 6000)



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

#####quick and dirty ilex
ilexflo<-filter(predy,GEN.SPA %in% c("ILE.MUC", "ILE.VER"))
ilexflo<-dplyr::select(ilexflo,GEN.SPA,Estimate,scenario)
colnames(ilexflo)[2]<-"Est.flo"

ilexbb<-filter(predy2,GEN.SPA %in% c("ILE.MUC", "ILE.VER"))
ilexbb<-dplyr::select(ilexbb,GEN.SPA,Estimate,scenario)
colnames(ilexbb)[2]<-"Est.bb"
ilex<-left_join(ilexflo,ilexbb)  
ilex$FLS<-ilex$Est.flo-ilex$Est.bb

setEPS()
postscript("Plots/climpredictions.eps",width = 12, height = 8)
predybig %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=phase,shape=phase),size=2.5)+geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0,color=phase),linetype="solid")+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=phase),linetype="dotted")+facet_wrap(~GEN.SPA,scale="free",ncol=3)+
  ggthemes::theme_base(base_size = 10)+scale_color_manual(values=c("darkgreen","darkorchid3"))+scale_shape_manual(values=c(17,19))
dev.off()


PHH<-data.frame[]

goo<-filter(predybig,GEN.SPA=="COR.COR")
goo%>%
arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=phase,shape=phase),size=2.5)+geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0,color=phase),linetype="solid",alpha=0.8)+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=phase),linetype="dotted",alpha=0.8)+
  ggthemes::theme_base(base_size = 10)+scale_color_manual(values=c("darkgreen","darkorchid3"))+scale_shape_manual(values=c(17,19))


jpeg("Plots/climpredictions_fixedscale.jpeg",width = 6, height = 10, units = 'in', res=300)
predybig %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=phase,shape=phase),size=2.5)+geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0,color=phase),linetype="solid",alpha=0.8)+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=phase),linetype="dotted",alpha=0.8)+facet_wrap(~GEN.SPA,scale="fixed",ncol=2)+theme_bw()+scale_color_manual(values=c("darkgreen","hotpink"))
dev.off()


HFchecker<-filter(HFreal, species %in%c('ACPE',"ACRU"))
ggplot(HFchecker)+geom_smooth(method="lm",aes(year,bb.jd),color="darkgreen",fill="darkgreen")+geom_smooth(method="lm",aes(year,fopn.jd),color="hotpink",fill="hotpink")+facet_wrap(~species)


modelhere <-mod.FLS.small
modelhere2<-mod.FLS.small
source("prep4plot.R")
source("exp_muplot_brms.R")
muplotfx(modelhere,modelhere2, "FLS change", 8, 8, c(0,7), c(-50, 120) , 130, 3.5)
dev.off()


### Vacccor is the only species that leafed and flower in every treatent


#### 
unique(dat$GEN.SPA)
hysters<-c("COM.PER","COR.COR","ACE.RUB")
hyst.dat<-filter(dat, GEN.SPA %in% hysters)

mod.hyst.flo<-brm(flo_day ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                   data = hyst.dat, family = gaussian(),control = list(adapt_delta = 0.99),
                   iter= 7000,
                   warmup = 6000)


ser.dat<-filter(dat, !GEN.SPA %in% hysters)
mod.ser.flo<-brm(flo_day ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                    data = ser.dat, family = gaussian(),control = list(adapt_delta = 0.99),
                    iter= 7000,
                    warmup = 6000)

spp <- unique(hyst.dat$GEN.SPA)
modelhere <-mod.hyst.flo
modelhere2<-mod.hyst.flo
source("prep4plot.R")
source("exp_muplot_brms.R")
muplotfx(modelhere,modelhere2, "Flo_hyst", 8, 8, c(0,7), c(-50, 120) , 130, 3.5)
dev.off()

spp <- unique(ser.dat$GEN.SPA)
modelhere <-mod.ser.flo
modelhere2<-mod.ser.flo
source("prep4plot.R")
source("exp_muplot_brms.R")
muplotfx(modelhere,modelhere2, "Flo_ser", 8, 8, c(0,7), c(-50, 120) , 130, 3.5)
dev.off()


save.image("new_flobud.mods.Rda")
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
