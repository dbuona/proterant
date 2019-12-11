rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_projs/")

library(ape)
library(phytools)
library(brms)
library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")
library("jpeg")

load("HFmodeloutput")
#useful function
extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}

##read in the data
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
HFsubber<-read.csv("HarvardForest/HFdata4modeling.csv",header=TRUE)
HF.tree<-read.tree("HarvardForest/HFtree4modeling.tre")
###make fls measure
HF$phys.fls<-HF$bb.jd-HF$fbb.jd
HF$funct.fls<-HF$l75.jd-HF$fopn.jd
HF$inter.fls<-HF$bb.jd-HF$fopn.jd

###make catagorical FLS
HF$hyst.funct<-ifelse(HF$funct.fls>0,1,0)
HF$hyst.phys<-ifelse(HF$phys.fls>0,1,0)
HF$hyst.inter<-ifelse(HF$inter.fls>0,1,0)

### prune the tree
HF<-dplyr::filter(HF,species!=("QUAL")) ## quercus alba has no flowers
spforcontmods<-HFsubber$species ##subset of species good for this analysis

HF.data<-dplyr::filter(HF,species %in% c(spforcontmods))
HF.data<-dplyr::left_join(HF.data,HFsubber, by="species") ###This is the data for the continuous models
##zscore predictors for these models
HF.data$pol_cent<-(HF.data$pol-mean(HF.data$pol,na.rm=TRUE))/(2*sd(HF.data$pol,na.rm=TRUE))
HF.data$precip_cent<-(HF.data$min_precip-mean(HF.data$min_precip))/(2*sd(HF.data$min_precip))
HF.data$flo_cent<-(HF.data$fopn.jd-mean(HF.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.data$fopn.jd,na.rm=TRUE))

HF.data$flo_cent.neg<--(HF.data$flo_cent)
HF.data$precip_cent.neg<--(HF.data$precip_cent)
###group by phylogeny
inv.phylo <- MCMCglmm::inverseA(HF.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###continuous models normal
modelcont.funct <- brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name), data = HF.data, 
                       family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 
pp_check(modelcont.funct)
###checks


#modelcont.phys <- brm(phys.fls~ pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent+(1|name), data = HF.data, 
 #                   family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 

#modelcont.inter<- brm(inter.fls~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
 #                     family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 



##binary
modelbin.funct<- brm(hyst.funct~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name), data = HF.data, 
                     family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000,warmup=3000) 

#modelbin.phys<- brm(hyst.phys~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent+(1|name), data = HF.data, 
 #                  family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000, warmup=3000)

#modelbin.inter<- brm(hyst.inter~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
       #             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000, warmup=3000)

funct.cont<-extract_coefs4HF(modelcont.funct)
funct.bin<-extract_coefs4HF(modelbin.funct)
funct.bin$class<-"functional"
funct.cont$class<-"functional"

#phys.cont<-extract_coefs4HF(modelcont.phys)
#phys.bin<-extract_coefs4HF(modelbin.phys)
#phys.cont$class<-"physiological"
#phys.bin$class<-"physiological"

#inter.cont<-extract_coefs4HF(modelcont.inter)
#inter.bin<-extract_coefs4HF(modelbin.inter)
#inter.cont$class<-"intermediate"
#inter.bin$class<-"intermediate"



cont<-funct.cont
bin<-funct.bin
cont$data_type<-"continuous"
bin$data_type<-"binary"


bin<-dplyr::filter(bin,trait!="Intercept")
bin$trait[which(bin$trait=="pol_cent")]<-"pollination syndrome"
bin$trait[which(bin$trait=="flo_cent")]<- "earlier flowering"
bin$trait[which(bin$trait=="precip_cent")]  <- "water dynamics"
bin$trait[which(bin$trait=="pol_cent:precip_cent")]<- "pollination:water dynamics"
bin$trait[which(bin$trait=="pol_cent:flo_cent")]<-"pollination:earlier flowering"
bin$trait[which(bin$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"
bin$data<-"HF"


cont<-dplyr::filter(cont,trait!="Intercept")
cont$trait[which(cont$trait=="pol")]<-"pollination syndrome"
cont$trait[which(cont$trait=="flo_cent")]<- "earlier flowering"
cont$trait[which(cont$trait=="precip_cent")]  <- "water dynamics"
cont$trait[which(cont$trait=="pol:precip_cent")]<- "pollination:water dynamics"
cont$trait[which(cont$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
cont$trait[which(cont$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"
cont$data<-"HF"


both<-rbind(cont,bin)
pd=position_dodgev(height=0.4)
b<-both %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(color=data),shape=1,position=pd,size=2,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=data),position=pd,width=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90, color=data),position=pd,width=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-30,30)+scale_color_manual(values=c("firebrick4"))+facet_wrap(~data_type)+
  labs(title = NULL,tag="b)")
 
#load(file = "MTSV_USFS/MTSVUSFS.mods")



##pd=position_dodgev(height=0.6)
#a<-comps %>%
#  arrange(estimate) %>%
#  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
#  ggplot(aes(estimate,trait))+geom_point(aes(shape=class,color=data,stroke=1.5),position=pd,size=2)+
#  geom_errorbarh(aes(xmin=low,xmax=high,linetype=class,color=data),position=pd,height=0)+
#  scale_linetype_manual(values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-8,8)+
#  scale_color_manual(values=c("orchid4","springgreen4"))+guides(size = "legend", linetype= "none")+
#  labs(title = NULL,tag="a)")

###sno var species or and bin
#meanhf<-HF.data %>% group_by(name,fopn.jd) %>% summarise(meanfunctFLS=mean(funct.fls,na.rm=TRUE))
#meanhf$FLSmeanfunctbin<-ifelse(meanhf$meanfunctFLS>0,1,0)
#meanhf<-left_join(meanhf,HFsubber)

#meanhf$pol_cent<-(meanhf$pol-mean(meanhf$pol,na.rm=TRUE))/(2*sd(meanhf$pol,na.rm=TRUE))
#meanhf$precip_cent<-(meanhf$min_precip-mean(meanhf$min_precip))/(2*sd(meanhf$min_precip))
#meanhf$flo_cent<-(meanhf$fopn.jd-mean(meanhf$fopn.jd,na.rm=TRUE))/(2*sd(meanhf$fopn.jd,na.rm=TRUE))

#meanhf$flo_cent.neg<--(meanhf$flo_cent)
#meanhf$precip_cent.neg<--(meanhf$precip_cent)

#modelbin.nosp<-brms::brm(FLSmeanfunctbin~ pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent+(1|name), data =meanhf, 
 #                    family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000,warmup=3000) 




(-.15)*(2*sd(HF.data$fbb.jd,na.rm=TRUE))+mean(HF.data$fbb.jd,na.rm=TRUE)

(max(HF.data$precip_cent)*(2*sd(HF.data$min_precip,na.rm=TRUE))+mean(HF.data$min_precip,na.rm=TRUE)) ##this is the maximum
(min(HF.data$precip_cent)*(2*sd(HF.data$min_precip,na.rm=TRUE))+mean(HF.data$min_precip,na.rm=TRUE)) ## this is minumom


(-.15)*(2*sd(HF.data$fbb.jd,na.rm=TRUE))+mean(HF.data$fbb.jd,na.rm=TRUE) ##Day 121 May first



#jpeg("muplots.jpeg",width = 7, height = 8, units = 'in', res=400)
#ggpubr::ggarrange(a+theme(axis.title=element_blank(),legend.title = element_blank() ),b+theme(axis.title.y=element_blank(),legend.title = element_blank() ),ncol=1,common.legend =FALSE, legend="top")
#dev.off()
apc.funct<- ggeffects::ggpredict(modelcont.funct,c("precip_cent","pol","flo_cent[-0.15]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot<-plot(apc.funct)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1,1.5),labels=c(6,13,19,26,33,40,47))+
  xlab("Min. precipitation across range (cm)")+ylab("Flowering to leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag="a)")+theme_linedraw()

#### unzscore precip values
1.5*(2*sd(HF.data$min_precip,na.rm=TRUE))+mean(HF.data$min_precip,na.rm=TRUE)

apc.funct.2<- ggeffects::ggpredict(modelcont.funct,c("flo_cent","pol","precip_cent[1]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot<-plot(apc.funct.2)+scale_x_continuous()+
  xlab("BOY flowering")+ylab("Flowering to leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag="a)")+theme_linedraw()




jpeg("HarvardForest/apcs.jpeg",width = 8.6, height = 4, units = 'in', res=200)
ggpubr::ggarrange(apc.funct,apc.funct.bin,common.legend = TRUE)
dev.off()


#Noww is Annual precipitation a good predictor
HF.weather<-read.csv("HarvardForest/mean.HF.precip.csv")
colnames(HF.weather)[1]<-"year"
HF.data<-left_join(HF.data,HF.weather, by="year")

HF.data$AP_cent<-(HF.data$AP-mean(HF.data$AP,na.rm=TRUE))/(2*sd(HF.data$AP,na.rm=TRUE))
HF.data$AP_cent.neg<--(HF.data$AP-mean(HF.data$AP,na.rm=TRUE))/(2*sd(HF.data$AP,na.rm=TRUE))


### Hysteranthous species
HFwind<-filter(HF.data,pol==1)


(1.5*(2*sd(HF.data$AP,na.rm=TRUE))+mean(HF.data$AP,na.rm=TRUE))

#fl;owering time or leafing time
summary(lm(funct.fls~fopn.jd,data=HF.data))
summary(lm(funct.fls~l75.jd,data=HF.data))

summary(lm(phys.fls~fbb.jd,data=HF.data))
summary(lm(phys.fls~bb.jd,data=HF.data))

##In hysteranthous speces
goo<-brm(funct.fls~precip_cent*AP_cent, data = HFwind, iter=4000, warmup=3000)
goo
save.image("HFmodeloutput")

