###HF modeling redo

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
library("phylolm")
library(ggstance)
load("finalmanuscriptmodels")

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

HFtree<-read.tree("HarvardForest/HFtree4modeling.tre" )

####mean only model
HF.tree<-drop.tip(HF.tree,tip = "Viburnum_lantanoides")
HF.tree$tip.label

HF.means<-HF.data %>% group_by(name)%>% summarise(meanFLS.func=mean(funct.fls,na.rm=TRUE),meanFLS.phys=mean(phys.fls,na.rm=TRUE),meanFLS.inter=mean(inter.fls,na.rm=TRUE),meanflo=mean(flo_cent,na.rm=TRUE),meandrought=mean(precip_cent),meanpol=mean(pol_cent))
HF.means$name

HF.means$hyster.fuct<-ifelse(HF.means$meanFLS.func>1,1,0)
HF.means$hyster.phys<-ifelse(HF.means$meanFLS.phys>1,1,0)
HF.means$hyster.inter<-ifelse(HF.means$meanFLS.inter>1,1,0)
HF.means<-  HF.means %>% remove_rownames %>% column_to_rownames(var="name")

z.phys.HFmean<-phylolm(meanFLS.phys~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean)

z.phys.HFbin<-phylolm(hyster.phys~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

z.funct.HFmean<-phylolm(meanFLS.func~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean)

z.funct.HFbin<-phylolm(hyster.fuct~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)


z.inter.HFmean<-phylolm(meanFLS.inter~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.inter.HFmean)

z.inter.HFbin<-phylolm(hyster.inter~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

inv.phylo <- MCMCglmm::inverseA(HF.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)




modelcont.funct.wspecies.ind<-brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                  family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 

modelcont.phys.wspecies.ind<-brm(phys.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

modelcont.inter.wspecies.ind<-brm(inter.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)


summary(z.funct.HFbin)
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
}
extract_CIs<-function(x){
  dplyr::filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
}
extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}

phys.cont<-extract_coefs4HF(modelcont.phys.wspecies.ind)
funct.cont<-extract_coefs4HF(modelcont.funct.wspecies.ind)

funct.cont$class<-"fopn-L75"
phys.cont$class<-"fbb-lbb"
cont<-rbind(phys.cont,funct.cont)

cont<-dplyr::filter(cont,trait!="Intercept")
cont$trait[which(cont$trait=="pol")]<-"pollination syndrome"
cont$trait[which(cont$trait=="flo_cent")]<- "earlier flowering"
cont$trait[which(cont$trait=="precip_cent")]  <- "water dynamics"
cont$trait[which(cont$trait=="pol:precip_cent")]<- "pollination:water dynamics"
cont$trait[which(cont$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
cont$trait[which(cont$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"


HF.funct.quan<-full_join(extract_coefs(z.funct.HFmean),extract_CIs(z.funct.HFmean),by="trait")

colnames(HF.funct.quan)<-c("trait","estimate","low","high")
HF.funct.quan$class<-"functional"

HF.phys.quan<-full_join(extract_coefs(z.phys.HFmean),extract_CIs(z.phys.HFmean),by="trait")
colnames(HF.phys.quan)<-c("trait","estimate","low","high")
HF.phys.quan$class<-"physiological"

quanty<-rbind(HF.phys.quan,HF.funct.quan)
quanty$model<-"quantitative" 

HF.funct.bin<-full_join(extract_coefs(z.funct.HFbin),extract_CIs(z.funct.HFbin),by="trait")

colnames(HF.funct.bin)<-c("trait","estimate","low","high")
HF.funct.bin$class<-"functional"

HF.phys.bin<-full_join(extract_coefs(z.phys.HFbin),extract_CIs(z.phys.HFbin),by="trait")
colnames(HF.phys.bin)<-c("trait","estimate","low","high")
HF.phys.bin$class<-"physiological"

binny<-rbind(HF.phys.bin,HF.funct.bin)
binny$model<-"categorical" 

comps<-rbind(binny,quanty)
comps<-dplyr::filter(comps,trait!="(Intercept)")
comps<-dplyr::filter(comps,trait!="sigma2")
unique(comps$trait)
comps$trait[which(comps$trait=="meanpol")]<-"pollination syndrome"
comps$trait[which(comps$trait=="meanflo")]<- "early flowering"
comps$trait[which(comps$trait=="meandrought")]  <- "water dynamics"
comps$trait[which(comps$trait=="meandrought:meanpol")]<- "pollination:water dynamics"
comps$trait[which(comps$trait=="meanflo:meanpol")]<-"pollination:early flowering"
comps$trait[which(comps$trait=="meanflo:meandrought")]<-"early flowering:water dynamics"
pd=position_dodgev(height=0.4)

physcomps<-filter(comps,class=="physiological")
functscomps<-filter(comps,class=="functional")


physcomps$class<-"fbb-lbb"
functscomps$class<-"fopn-L75"
comps<-rbind(physcomps,functscomps)

binary.main<-filter(physcomps,model=="categorical")
quant.main<-filter(physcomps,model!="categorical")

a<-binary.main %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(),position=pd,size=2,color="firebrick1")+
  geom_errorbarh(aes(xmin=low,xmax=high),position=pd,height=0,size=.5,color="firebrick1")+
  xlim(-1.5,2)+ theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
 guides(size = "legend", linetype= "none")+theme(axis.title.x=element_blank(),
                                                                plot.margin = margin(l=5,r = 1,t=20,b=10))


b<-quant.main %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(),position=pd,size=2,color="goldenrod3")+
  geom_errorbarh(aes(xmin=low,xmax=high),position=pd,height=0,size=.5,color="goldenrod3")+
  xlim(-35,50)+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
 guides(size = "legend", linetype= "none")+theme(axis.title.y=element_blank(),
                                                                axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.x=element_blank(),
                                                                plot.margin = margin(r = 1.5, l = 1.5,t=20,b=10))
                                                                

cont.phys<-filter(cont,class=="fbb-lbb")
c<-cont.phys %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(),position=pd,size=2,stroke=.5,color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,height=0,linetype="dotted",color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,height=0,linetype="solid",color="darkorchid3")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-35,50)+theme(axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.title.x=element_blank(),
                     axis.ticks.y=element_blank(),
                                    plot.margin = margin(l = 1,r=5,t=20,b=10) )



ggpubr::ggarrange(a,b,c,nrow=1,ncol=3,widths=c(2,1,1),labels=c("a)","b)","c)"),hjust=c(-11.2,-.5,-.5))
?ggarrange()



######Suppliment

bin.comps<-filter(comps,model=="categorical")
quant.comps<-filter(comps,model!="categorical")

supa<-bin.comps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class),position=pd,size=2,color="firebrick1")+
  geom_errorbarh(aes(xmin=low,xmax=high,group=class),position=pd,height=0,size=.5,color="firebrick1")+
  xlim(-1.5,2)+ theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  guides(size = "legend", linetype= "none")+theme(axis.title.x=element_blank(),
                                                  plot.margin = margin(l=5,r = 1,t=20,b=10))

supb<-quant.comps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class),position=pd,size=2,color="goldenrod3")+
  geom_errorbarh(aes(xmin=low,xmax=high,group=class),position=pd,height=0,size=.5,color="goldenrod3")+
  xlim(-45,50)+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  guides(size = "legend", linetype= "none")+theme(axis.title.y=element_blank(),
                                                  axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.x=element_blank(),
                                                  plot.margin = margin(r = 1.5, l = 1.5,t=20,b=10))

supc<-cont %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=class),position=pd,size=2,stroke=.5,color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,height=0,linetype="dotted",color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,height=0,linetype="solid",color="darkorchid3")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-45,50)+theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.y=element_blank(),
                     plot.margin = margin(l = 1,r=5,t=20,b=10) )

ggpubr::ggarrange(supa,supb,supc,nrow=1,ncol=3,widths=c(2,1,1),labels=c("a)","b)","c)"),hjust=c(-11.2,-.5,-.5),common.legend = TRUE,legend = "right")

data.frame(phenophase=c("leaf budburst","flower budburst","flowers open", "75% of leaves full size"),HF=c("lbb","fbb","fopn","l75"),BBCH=c("09","55","60","17"))

##### alternative models
colnames(HF.data)


modelcont.phys.wspecies.ind<-brm(phys.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)
save.image("HFfinalmanuscriptmods")

###### alternative models:
HF.data$dummydrought=NA
HF.data$dummydrought[which(HF.data$drought_tol=="H")]<-1
HF.data$dummydrought[which(HF.data$drought_tol=="M")]<-.5
HF.data$dummydrought[which(HF.data$drought_tol=="L")]<-0
HF.data$dummy_cent<-(HF.data$dummydrought-mean(HF.data$dummydrought,na.rm=TRUE))/(2*sd(HF.data$dummydrought,na.rm=TRUE))


modelcont.phys.dummydrought<-brm(phys.fls~ pol+flo_cent+dummy_cent+dummy_cent:flo_cent+dummy_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

save.image("finalmanuscriptmodels") 
