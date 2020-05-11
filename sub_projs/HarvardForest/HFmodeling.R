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
#modelcont.funct <- brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name), data = HF.data, 
 #                      family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 
###checks

modelcont.funct.wspecies.ind<-brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 

modelcont.phys.wspecies.ind<-brm(phys.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                  family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 
pp_check(modelcont.phys.wspecies.ind)


prior <- c(prior(student_t(3, 16, 45), class = Intercept),
           prior(student_t(3, 0, 45) , class = sd),
           prior(student_t(3, 0, 45) , class = sigma))

modelcont.funct.wspecies.ind.altprior<-brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                  family = gaussian(), cov_ranef = list(name= A),prior=prior,control=list(adapt_delta=0.95),iter=4000, warmup=3000) 




##binary
modelbin.funct.sps.ind<- brm(hyst.funct~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                     family = bernoulli(link="logit"), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000,warmup=3000) 

modelbin.funct.sps<- brm(hyst.funct~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|species), data = HF.data, 
                             family = bernoulli(link="logit"), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000,warmup=3000) 

modelbin.phys.sps.ind<- brm(hyst.phys~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                             family = bernoulli(link="logit"), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000,warmup=3000) 

prior_summary(modelcont.funct.wspecies.ind)
prior_summary(modelbin.funct.sps)

tau2 <- brms::VarCorr(modelbin.funct.wspecies.ind)[[1]]$sd[1]^2



#####phylo signal
# computing the ICC for the intercept
ICC1 <- tau2 / (tau2 + (pi^2 / 3) )
### phylo signal
hyp <- "sd_name__Intercept^2 / (sd_name__Intercept^2 + sigma^2) = 0"

lambda.FLS <- hypothesis(modelcont.funct.wspecies.ind, hyp, class = NULL)

par(mfrow=c(1,2))
d2 <- density(lambda.FLS$samples[,1])
plot(d2,main="",xlab="",
     xlim=c(0,1),col="darkblue")
polygon(d2, col=adjustcolor("darkblue",0.4), border="darkblue")
abline(v=mean(lambda.FLS$samples[,1]),lty=2,col="blue")


hyp2<-"sd_name__Intercept^2/(sd_name__Intercept^2+ ( 3.141593^2 / 3)) = 0"

lambdabin<-hypothesis(modelbin.funct, hyp2, class = NULL)

d <- density(lambdabin$samples[,1])
plot(d,main="",xlab="",
     xlim=c(0,1),col="darkblue")
polygon(d, col=adjustcolor("darkblue",0.4), border="darkblue")
abline(v=mean(lambdabin$samples[,1]),lty=2,col="blue")



###alternative phenophases
#modelbin.phys<- brm(hyst.phys~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent+(1|name), data = HF.data, 
 #                  family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000, warmup=3000)

#modelbin.inter<- brm(hyst.inter~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.data, 
       #             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000, warmup=3000)

funct.cont<-extract_coefs4HF(modelcont.funct.wspecies.ind)
funct.bin<-extract_coefs4HF(modelbin.funct.sps.ind)

phys.cont<-extract_coefs4HF(modelcont.phys.wspecies.ind)
phys.bin<-extract_coefs4HF(modelbin.phys.sps.ind)
phys.cont$class<-"fbb-lbb"
funct.cont$class="fopn-L75"
phys.bin$class<-"fbb-lbb"
funct.bin$class<-"fopn-L75"
#inter.cont<-extract_coefs4HF(modelcont.inter)
#inter.bin<-extract_coefs4HF(modelbin.inter)
#inter.cont$class<-"intermediate"
#inter.bin$class<-"intermediate"



cont<-rbind(funct.cont,phys.cont)
bin<-rbind(funct.bin,phys.bin)
cont$data_type<-"continuous"
bin$data_type<-"categorical"


bin<-dplyr::filter(bin,trait!="Intercept")
bin$trait[which(bin$trait=="pol")]<-"pollination syndrome"
bin$trait[which(bin$trait=="flo_cent")]<- "earlier flowering"
bin$trait[which(bin$trait=="precip_cent")]  <- "water dynamics"
bin$trait[which(bin$trait=="pol:precip_cent")]<- "pollination:water dynamics"
bin$trait[which(bin$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
bin$trait[which(bin$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"



cont<-dplyr::filter(cont,trait!="Intercept")
cont$trait[which(cont$trait=="pol")]<-"pollination syndrome"
cont$trait[which(cont$trait=="flo_cent")]<- "earlier flowering"
cont$trait[which(cont$trait=="precip_cent")]  <- "water dynamics"
cont$trait[which(cont$trait=="pol:precip_cent")]<- "pollination:water dynamics"
cont$trait[which(cont$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
cont$trait[which(cont$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"

cont %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=class),position=pd,size=3,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,height=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-40,30)+scale_color_manual(values=c("firebrick4"))+scale_shape_discrete(name = "data type")

bin %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=class),position=pd,size=3,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,height=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-15,10)+scale_color_manual(values=c("firebrick4"))+scale_shape_discrete(name = "data type")



both<-rbind(cont,bin)
#jpeg("HF.jpeg",width = 6, height = 4, units = 'in', res=350)
tiff("HF.tiff",width = 6, height = 4, units = 'in', res=350)
pd=position_dodgev(height=0.4)
both %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=data_type),position=pd,size=3,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=data_type),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=data_type),position=pd,height=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-30,30)+scale_color_manual(values=c("firebrick4"))+scale_shape_discrete(name = "data type")
dev.off() 
load(file = "MTSV_USFS/MTSVUSFS.mods")



pd=position_dodgev(height=0.6)
a<-comps %>%
  arrange(estimate) %>%
 mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class,color=data,stroke=1.5),position=pd,size=2)+
  geom_errorbarh(aes(xmin=low,xmax=high,linetype=class,color=data),position=pd,height=0,size=1)+
  scale_linetype_manual(values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-8,9)+
  scale_color_manual(values=c("orchid4","springgreen4"))+guides(size = "legend", linetype= "none")+
  labs(title = NULL,tag="a)")

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


(.51)*(2*sd(HF.data$fbb.jd,na.rm=TRUE))+mean(HF.data$fbb.jd,na.rm=TRUE) ##Day 121 May first



jpeg("muplots.jpeg",width = 7, height = 8, units = 'in', res=400)
ggpubr::ggarrange(a+theme(axis.title=element_blank(),legend.title = element_blank() ),b+theme(axis.title.y=element_blank(),legend.title = element_blank() ),widths=c(1,2),ncol=2,common.legend =FALSE, legend="top")
dev.off()

apc.funct<- ggeffects::ggpredict(modelcont.funct.wspecies.ind,c("precip_cent","pol","flo_cent[-0.15]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot<-plot(apc.funct)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1,1.5),labels=c(6,13,19,26,33,40,47))+
  xlab("Min. precipitation across range (cm)")+ylab("Flowering to leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag=NULL)+theme_linedraw()

apc.funct2<- ggeffects::ggpredict(modelcont.funct.wspecies.ind,c("precip_cent","pol","flo_cent[.55]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot2<-plot(apc.funct2)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1,1.5),labels=c(6,13,19,26,33,40,47))+
  xlab("Min. precipitation across range (cm)")+ylab("Flowering to leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag=NULL)+theme_linedraw()





#jpeg("HarvardForest/apcs.jpeg",width = 8.6, height = 4, units = 'in', res=200)
tiff("HarvardForest/apcs.tiff",width = 6, height = 4, units = 'in', res=300)
apc.funct.plot
dev.off()


#Noww is Annual precipitation a good predictor
HF.weather<-read.csv("HarvardForest/mean.HF.precip.csv")
colnames(HF.weather)[1]<-"year"
HF.data<-left_join(HF.data,HF.weather, by="year")

HF.data$AP_cent<-(HF.data$AP-mean(HF.data$AP,na.rm=TRUE))/(2*sd(HF.data$AP,na.rm=TRUE))

AP.mod<-brm(funct.fls~ AP, data = HF.data, family = gaussian() ,iter=5000, warmup=3000) 
summary(AP.mod)



goo2<-brm(funct.fls~ AP*pol, data = HF.data, family = gaussian() ,iter=3000, warmup=2000) 
summary(goo2)



##### alternative FLS models

modelcont.phys.wspecies.ind<-brm(phys.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                  family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 

modelcont.inter.wspecies.ind<-brm(inter.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 

modelbin.phys.sps<- brm(hyst.phys~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                             family = bernoulli(link="logit"), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000,warmup=3000) 

modelbin.inter.sps.ind<- brm(hyst.inter~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                            family = bernoulli(link="logit"), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000,warmup=3000) 



funct.cont<-extract_coefs4HF(modelcont.phys.wspecies.ind)
funct.bin<-extract_coefs4HF(modelbin.phys.sps.ind)


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
bin$data_type<-"categorical"


bin<-dplyr::filter(bin,trait!="Intercept")
bin$trait[which(bin$trait=="pol")]<-"pollination syndrome"
bin$trait[which(bin$trait=="flo_cent")]<- "earlier flowering"
bin$trait[which(bin$trait=="precip_cent")]  <- "water dynamics"
bin$trait[which(bin$trait=="pol:precip_cent")]<- "pollination:water dynamics"
bin$trait[which(bin$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
bin$trait[which(bin$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"



cont<-dplyr::filter(cont,trait!="Intercept")
cont$trait[which(cont$trait=="pol")]<-"pollination syndrome"
cont$trait[which(cont$trait=="flo_cent")]<- "earlier flowering"
cont$trait[which(cont$trait=="precip_cent")]  <- "water dynamics"
cont$trait[which(cont$trait=="pol:precip_cent")]<- "pollination:water dynamics"
cont$trait[which(cont$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
cont$trait[which(cont$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"




both<-rbind(cont,bin)
jpeg("HF.phys.jpeg",width = 6, height = 4, units = 'in', res=350)
pd=position_dodgev(height=0.4)
both %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=data_type),position=pd,size=3,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=data_type),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=data_type),position=pd,height=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-35,30)+scale_color_manual(values=c("firebrick4"))+scale_shape_discrete(name = "data type")
dev.off() 


load(file = "MTSV_USFS/MTSVUSFS.mods")

save.image("HFmodeloutput")


