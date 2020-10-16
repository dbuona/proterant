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

dat<-read.csv("flobudsdata.use.csv",header = TRUE)

dat$Light<-ifelse(dat$Light=="S",0,1)
dat$Force<-ifelse(dat$Force=="C",0,1)

dat<-filter(dat, !GEN.SPA %in% c("AME.SPP","BET.SPP", "BET.ALL","ACE.SAC"))

###filter data
dat.chillmet<-filter(dat,Chill==1)
dat.chillmet<-filter(dat.chillmet,Light==1)

#####plasticity models for 3 bbch stages 


mod.bb.int.phh<-brm(budburst.9.~Force+(Force|GEN.SPA),
                data = dat.chillmet, family = gaussian(),
                iter= 4000,
                warmup = 3000)  




mod.flo.int<-brm(flo_day~Force+(Force|GEN.SPA),
                 data = dat.chillmet, family = gaussian(),
                 iter= 4000,
                 warmup = 3000) 


coef(mod.flo.int)
coef(mod.bb.int.phh)

extract_coefs<-function(x){rownames_to_column(as.data.frame(coef(x, summary=TRUE,probs=c(0.25,0.75))),"GEN.SPA")
}


flo<-extract_coefs(mod.flo.int)
flo$phase<-"reproductive"
leaf<-extract_coefs(mod.bb.int.phh)
leaf$phase<-"vegetative"
phh<-rbind(flo,leaf)
colnames(phh)

phh$order<-NA

phh$order[which(phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase=="reproductive")]<-"first"
phh$order[which(phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase!="reproductive")]<-"second"

phh$order[which(!phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase=="reproductive")]<-"second"

phh$order[which(!phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase!="reproductive")]<-"first"
phh2<-dplyr::select(phh,GEN.SPA,GEN.SPA.Estimate.Force,GEN.SPA.Est.Error.Force,GEN.SPA.Q25.Force,GEN.SPA.Q75.Force,phase,order)
colnames(phh2)<-c("Species","Estimate","error","Q25","Q75","phase","sequence")
phh2<-phh2 %>%arrange(Species,sequence)

save.image("phh.mod.output.Rda")
png("Plots/Flobuds_manuscript_figs/phh_plot.png",width = 5,height = 5,units = "in",res=300)
ggplot(phh,aes(GEN.SPA.Estimate.Force,GEN.SPA))+geom_point(aes(shape=phase,color=order),size=2)+
geom_errorbarh(aes(xmin=GEN.SPA.Q25.Force,xmax=GEN.SPA.Q75.Force,group=phase,color=order),height=0)+scale_color_brewer(type="qual",palette = 2)+  
geom_vline(xintercept = 0,linetype="dashed")+ggthemes::theme_base(base_size = 11)+ylab("Species")+xlab("Sensitivity to forcing")
dev.off()
