###TO DO run models without ACE.SAC and BET.ALL to see if figures work better
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
load("writing/phh.mod.output.Rda")
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
load("phh.mod.output.Rda")
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

extract_coefs<-function(x){rownames_to_column(as.data.frame(coef(x, summary=TRUE,probs=c(0.025,0.25,0.75,0.975))),"GEN.SPA")
}


flo<-extract_coefs(mod.flo.int)
flo$phase<-"flowering"
leaf<-extract_coefs(mod.bb.int.phh)
leaf$phase<-"leaf budburst"
phh<-rbind(flo,leaf)
colnames(phh)

phh$order<-NA

phh$order[which(phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase=="flowering")]<-"first"
phh$order[which(phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase!="flowering")]<-"second"

phh$order[which(!phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase!="leaf budburst")]<-"second"

phh$order[which(!phh$GEN.SPA %in% c("COM.PER","ACE.RUB","COR.COR") & phh$phase=="leaf budburst")]<-"first"

phh2<-dplyr::select(phh,GEN.SPA,GEN.SPA.Estimate.Force,GEN.SPA.Est.Error.Force,GEN.SPA.Q2.5.Force,GEN.SPA.Q25.Force,GEN.SPA.Q75.Force,GEN.SPA.Q97.5.Force,phase,order)
colnames(phh2)<-c("Species","Estimate","error","Q2.5","Q25","Q75","Q97.5","phase","sequence")
phh2<-phh2 %>%arrange(Species,sequence)

GEN.SPA<-unique(phh$GEN.SPA)
Species<-c("Acer pensylvanicum","Acer rubrum","Comptonia perigrina",
           "Corylus cornuta","Ilex mucronata", "Ilex verticillata", "Prunus penylvanica","Prunus virginiana",
           "Vaccinium corymbosum","Viburnum acerifolium")
coresp<-data.frame(GEN.SPA,Species)

phh<-left_join(phh,coresp)
save.image("phh.mod.output.Rda")
png("Plots/Flobuds_manuscript_figs/phh_plot.png",width = 5,height = 3,units = "in",res=250)
ggplot(phh,aes(GEN.SPA.Estimate.Force,Species))+geom_point(aes(shape=phase,color=order),size=2)+
geom_errorbarh(aes(xmin=GEN.SPA.Q25.Force,xmax=GEN.SPA.Q75.Force,group=phase,color=order),height=0)+scale_color_brewer(type="qual",palette = 2)+  
geom_vline(xintercept = 0,linetype="dashed")+ggthemes::theme_few(base_size = 11)+ylab("Species")+xlab("Sensitivity to forcing")+
  theme(axis.text.y = element_text(face="italic"))
dev.off()
