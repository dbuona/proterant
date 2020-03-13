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

load("uhohmods")
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


##### models with mean flower time
#a) with and without phylogeny
meanflo<-HF.data %>% group_by(name) %>% summarise(meanflotime=mean(fopn.jd,na.rm=TRUE))
HF.data<-left_join(HF.data,meanflo)

meanFLS<-HF.data %>% group_by(name) %>% summarise(meanFLS=mean(funct.fls,na.rm=TRUE))
HF.data<-left_join(HF.data,meanFLS)

HF.data$meanflocent<-(HF.data$meanflotime-mean(HF.data$meanflotime,na.rm=TRUE))/(2*sd(HF.data$meanflotime,na.rm=TRUE))

##no variationa t all model
modelcont.phylo.means<-brm(meanFLS~ pol+meanflocent+precip_cent+(1|name), data = HF.data, 
                         family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99),iter=8000, warmup=7000) ### same

+precip_cent:meanflocent+precip_cent:pol+pol:meanflocent+

###models par A
modelcont.nophylo.ind<-brm(funct.fls~ pol+meanflocent+precip_cent+precip_cent:meanflocent+precip_cent:pol+pol:meanflocent+(1|species), data = HF.data, 
                         family = gaussian(),control=list(adapt_delta=0.99),iter=8000, warmup=7000) ##ESS on sd is very low

modelcont.phylo.ind<-brm(funct.fls~ pol+meanflocent+precip_cent+precip_cent:meanflocent+precip_cent:pol+pol:meanflocent+(1|name)+(1|species), data = HF.data, 
                                         family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99),iter=8000, warmup=7000) ### same

##Full with phylo
HF.data$varflo <- HF.data$fopn.jd-HF.data$meanflotime 
HF.data$varflo_cent<-(HF.data$varflo-mean(HF.data$varflo,na.rm=TRUE))/(2*sd(HF.data$varflo,na.rm=TRUE))

modelcont.var.sp.phylo<-brm(funct.fls~ pol+meanflocent+varflo_cent+precip_cent+precip_cent:meanflocent+precip_cent:pol+pol:meanflocent+(1|name)+(1|species), data = HF.data, 
                         family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99),iter=8000, warmup=7000) ### same

modelcont.var.sp.phylo2<-brm(inter.fls~ pol+meanflocent+varflo_cent+precip_cent+precip_cent:meanflocent+precip_cent:pol+pol:meanflocent+(1|name)+(1|species), data = HF.data, 
                            family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99),iter=8000, warmup=7000)

modelphys.var.sp.phylo2<-brm(phys.fls~ pol+meanflocent+varflo_cent+precip_cent+precip_cent:meanflocent+precip_cent:pol+pol:meanflocent+(1|name)+(1|species), data = HF.data, 
                             family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99),iter=8000, warmup=7000)


apc.funct<- ggeffects::ggpredict(modelcont.var.sp.phylo,c("precip_cent","pol","varflo_cent[0]","meanflocent[-1.5,.5]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot<-plot(apc.funct)
+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1,1.5),labels=c(6,13,19,26,33,40,47))+
  xlab("Min. precipitation across range (cm)")+ylab("Flowering to leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag=NULL)+theme_linedraw()


mean.phylo<-extract_coefs4HF(modelcont.phylo.ind)
mean.no.phylo<-extract_coefs4HF(modelcont.nophylo.ind)


mean.phylo$data_type<-"phylo"
mean.no.phylo$data_type<-"no phylo"

cont<-rbind(mean.no.phylo,mean.phylo)
cont$approach<-"novar predictors"

cont<-dplyr::filter(cont,trait!="Intercept")
cont$trait[which(cont$trait=="pol")]<-"pollination syndrome"
cont$trait[which(cont$trait=="meanflocent")]<- "earlier flowering"
cont$trait[which(cont$trait=="precip_cent")]  <- "water dynamics"
cont$trait[which(cont$trait=="pol:precip_cent")]<- "pollination:water dynamics"
cont$trait[which(cont$trait=="pol:meanflocent")]<-"pollination:earlier flowering"
cont$trait[which(cont$trait=="meanflocent:precip_cent")]<-"earlier flowering:water dynamics"



##part b
modelcont.funct.var.nophylo<-brm(inter.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|tree.id/species),
                                  data = HF.data, family = gaussian(),control=list(adapt_delta=0.99),iter=7000, warmup=5000)



modelcont.funct.var.phylo<-brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                        family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 


var.phylo<-extract_coefs4HF(modelcont.funct.var.phylo)
var.no.phylo<-extract_coefs4HF(modelcont.funct.var.nophylo)


var.phylo$data_type<-"phylo"
var.no.phylo$data_type<-"no phylo"

cont2<-rbind(var.no.phylo,var.phylo)
cont2$approach<-"vary predictors"

cont2<-dplyr::filter(cont2,trait!="Intercept")
cont2$trait[which(cont2$trait=="pol")]<-"pollination syndrome"
cont2$trait[which(cont2$trait=="flo_cent")]<- "earlier flowering"
cont2$trait[which(cont2$trait=="precip_cent")]  <- "water dynamics"
cont2$trait[which(cont2$trait=="pol:precip_cent")]<- "pollination:water dynamics"
cont2$trait[which(cont2$trait=="pol:flo_cent")]<-"pollination:earlier flowering"
cont2$trait[which(cont2$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"

both<-rbind(cont,cont2)
both<-tidyr::unite(both,merger,data_type,approach,remove=FALSE)
?unite()
pd=position_dodgev(height=0.4)
jpeg("Uhoh_best.jpeg",width = 6, height = 3, units = 'in', res=300)
both %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","pollination:variation flowering","earlier flowering","variation: flowering","water dynamics","pollination syndrome"))) %>% 
  ggplot(aes(Estimate,trait))+geom_point(aes(color=approach,shape=data_type,group=merger),position=pd,size=3,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=approach,group=merger),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=approach,group=merger),position=pd,height=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+scale_shape_discrete(name = "data type")
dev.off()


save.image("uhohmods")




