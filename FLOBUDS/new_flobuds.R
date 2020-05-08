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

###plot raw dat. nor sure these plot are so helpful. maybe better to do predictions
#plot.raw.dat<-gather(dat,phase,DOY,2:4)
#plot.raw.dat<-unite(plot.raw.dat,Treatcode,5:7,remove=FALSE)
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_1_0"]<-"WL28"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_1_1"]<-"WL56"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_1_0"]<-"CL28"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_0_0"]<-"CS28"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_0_1"]<-"CS56"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_1_1"]<-"CL56"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_0_1"]<-"WS56"
#plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_0_0"]<-"WS28"
#plot.raw.dat$phase[plot.raw.dat$phase=="budburst.9."]<-"budburst"
#plot.raw.dat$phase[plot.raw.dat$phase=="leaf_day.15."]<-"expansion"
#plot.raw.dat$phase[plot.raw.dat$phase=="flo_day"]<-"flower"

#plot.raw.dat.lo<-filter(plot.raw.dat,phase!="budburst")
#plot.raw.dat.bb<-filter(plot.raw.dat,phase!="expansion")
#goodsp<-c("ACE.PEN","COM.PER","COR.COR","ILE.MUC","PRU.PEN","VAC.COR")
#plot.raw.dat.lo<-filter(plot.raw.dat.lo,GEN.SPA %in% c(goodsp) )
#plot.raw.dat.bb<-filter(plot.raw.dat.bb,GEN.SPA %in% c(goodsp) )
##clean up variables

#plot.raw.dat$taxa<-NA
#plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="ACE.PEN"]<-"A. pensylvanicum"
#plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="COM.PER"]<-"C. peregrina"
#plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="COR.COR"]<-"C. cornuta"
#plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="ILE.MUC"]<-"I. mucronata"
#plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="PRU.PEN"]<-"P. pensylvanica"
#plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="VAC.COR"]<-"V corymbosum"

#goodspfull<-filter(plot.raw.dat,GEN.SPA %in% c(goodsp)) 
#goodspfull$force<-ifelse(goodspfull$Force==1,"High","Low")
#goodspfull$photo<-ifelse(goodspfull$Light==1,"Long day","Short day")
#goodspfull$chill<-ifelse(goodspfull$Chill==1,"High chill","Low chill")

#ggplot(goodspfull,aes(taxa,DOY))+geom_point(aes(color=phase,shape=force),size=0.8)+stat_summary(aes(color=phase,shape=force))+facet_grid(chill~photo)+theme_base()+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Treatment combination")+theme_base(base_size = 7)

#pdf("Plots/flo_buds_figures/goodsps_rawplots2.pdf",width=11, height=6)
#ggplot(goodspfull,aes(Treatcode,DOY))+geom_point(aes(color=phase,shape=phase),size=.8)+stat_summary(aes(color=phase,shape=phase))+facet_wrap(~taxa)+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Treatment combination")+theme_base(base_size = 7)
#pd=position_dodge2(width=.3,preserve="single")
#?position_dodge2()
#ggplot(goodspfull,aes(force,DOY))+geom_point(aes(shape=photo,color=phase),size=0.8)+stat_summary(aes(shape=photo,color=phase),position=pd)+facet_grid(chill~taxa)+theme_base()+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Forcing")+theme_base(base_size = 12)
#dev.off()


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



muplotfx(modelhere,modelhere2, "budburst vs. flowering", 8, 8, c(0,8), c(-40, 30) , 40, 3.5,40,4.5)
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



setEPS()
postscript("Plots/climpredictions.eps",width = 12, height = 8)
predybig %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=phase,shape=phase),size=2.5)+geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0,color=phase),linetype="solid")+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=phase),linetype="dotted")+facet_wrap(~GEN.SPA,scale="free",ncol=5)+
  ggthemes::theme_base(base_size = 10)+scale_color_manual(values=c("darkgreen","darkorchid3"))+scale_shape_manual(values=c(17,19))
dev.off()

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
