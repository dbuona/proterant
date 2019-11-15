rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")

library(ape)
library(phytools)
library(brms)
library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")

#useful function
extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}

##read in the data
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)

###make fls measure
HF$phys.fls<-HF$bb.jd-HF$fbb.jd
HF$funct.fls<-HF$l75.jd-HF$fopn.jd
HF$inter.fls<-HF$bb.jd-HF$fopn.jd

###make catagorical FLS
HF$hyst.funct<-ifelse(HF$funct.fls>0,1,0)
HF$hyst.phys<-ifelse(HF$phys.fls>0,1,0)
HF$hyst.inter<-ifelse(HF$inter.fls>0,1,0)

source("..//Scripts/continuous_mod_prep.R") ### prune the tree
HF<-dplyr::filter(HF,species!=("QUAL")) ## quercus alba has no flowers
spforcontmods<-df$species ##subset of species good for this analysis

HF.data<-dplyr::filter(HF,species %in% c(spforcontmods))
HF.data<-dplyr::left_join(HF.data,traits, by="species") ###This is the data for the continuous models
##zscore predictors for these models
HF.data$pol_cent<-(HF.data$pol-mean(HF.data$pol,na.rm=TRUE))/(2*sd(HF.data$pol,na.rm=TRUE))
HF.data$precip_cent<-(HF.data$min_precip-mean(HF.data$min_precip))/(2*sd(HF.data$min_precip))
HF.data$flo_cent<-(HF.data$fopn.jd-mean(HF.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.data$fopn.jd,na.rm=TRUE))
HF.data$flo_cent.neg<--(HF.data$flo_cent)

###group by phylogeny
inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###continuous models
modelcont.funct <- brm(funct.fls~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
                       family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 


modelcont.phys <- brm(phys.fls~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
                      family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 

modelcont.inter<- brm(inter.fls~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
                      family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 



##binary
modelbin.funct<- brm(hyst.funct~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
                     family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000,warmup=3000) 

modelbin.phys<- brm(hyst.phys~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
                    family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000, warmup=3000)

modelbin.inter<- brm(hyst.inter~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
                    family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000, warmup=3000)

funct.cont<-extract_coefs4HF(modelcont.funct)
funct.bin<-extract_coefs4HF(modelbin.funct)
funct.bin$class<-"functional"
funct.cont$class<-"functional"

phys.cont<-extract_coefs4HF(modelcont.phys)
phys.bin<-extract_coefs4HF(modelbin.phys)
phys.cont$class<-"physiological"
phys.bin$class<-"physiological"

inter.cont<-extract_coefs4HF(modelcont.inter)
inter.bin<-extract_coefs4HF(modelbin.inter)
inter.cont$class<-"intermediate"
inter.bin$class<-"intermediate"



cont<-inter.cont
bin<-inter.bin
cont$data_type<-"continuous"
bin$data_type<-"binary"


bin<-dplyr::filter(bin,trait!="Intercept")
bin$trait[which(bin$trait=="pol_cent")]<-"pollination syndrome"
bin$trait[which(bin$trait=="flo_cent.neg")]<- "earlier flowering"
bin$trait[which(bin$trait=="precip_cent")]  <- "water dynamics"
bin$trait[which(bin$trait=="pol_cent:precip_cent")]<- "pollination:water dynamics"
bin$trait[which(bin$trait=="pol_cent:flo_cent.neg")]<-"pollination:flowering"
bin$trait[which(bin$trait=="flo_cent.neg:precip_cent")]<-"flowering:water dynamics"


pd=position_dodgev(height=0.4)
binplot<-bin %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=class),position=pd,size=3,stroke=1.5)+scale_shape_manual(name="data type",values=c(21,22,23))+#scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,width=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,width=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")


cont<-dplyr::filter(cont,trait!="Intercept")
cont$trait[which(cont$trait=="pol_cent")]<-"pollination syndrome"
cont$trait[which(cont$trait=="flo_cent.neg")]<- "earlier flowering"
cont$trait[which(cont$trait=="precip_cent")]  <- "water dynamics"
cont$trait[which(cont$trait=="pol_cent:precip_cent")]<- "pollination:water dynamics"
cont$trait[which(cont$trait=="pol_cent:flo_cent.neg")]<-"pollination:flowering"
cont$trait[which(cont$trait=="flo_cent.neg:precip_cent")]<-"flowering:water dynamics"


contplot<-cont %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=class),position=pd,size=3,stroke=1.5)+scale_shape_manual(name="data type",values=c(21,22,23))+#scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,width=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,width=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")

ggpubr::ggarrange(binplot,contplot,common.legend = TRUE)

### no show only interactions with 80% CI not overlating zero 

##take a model where there are interactions

goober<- ggpredict(modelbin.phys,c("precip_cent","pol_cent","flo_cent.neg[.2]"), ci.lvl=0.50) ##May 15
apc.phys<-plot(goober)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1),labels=c(15,33,50,67,84,101))+
  xlab("Min. precipitation across range (cm)")+ylab("Likelihood of flower before leaf budburst")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag="a)")+theme_linedraw()

goober2<- ggpredict(modelcont.funct,c("precip_cent","pol_cent","flo_cent.neg[.2]"), ci.lvl=0.50)  #may 15
apc.funct<-plot(goober2)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1),labels=c(15,33,50,67,84,101))+
  xlab("Min. precipitation across range (cm)")+ylab("Flowers opening to 75% leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag="b)")+theme_linedraw()



