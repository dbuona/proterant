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
#load("finalmanuscriptmodels")
load("altpredictors") 

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


####begin the model
HFtree<-read.tree("HarvardForest/HFtree4modeling.tre" )

#### prune the tree to the same species
HF.tree<-drop.tip(HF.tree,tip = "Viburnum_lantanoides")
HF.tree$tip.label

###add aridity index:
arid<-read.csv("..//Data/arid_indexranges.csv")
arid$Taxon.name<-sub(" ", "_",arid$Taxon.name , fixed=TRUE)
colnames(arid)[1]<-"name"
arid$name<-sub(" ", "", arid$name)

intersect(unique(HF.data$name),arid$name)
##impute missing data
addos<-data.frame(name=setdiff(unique(HF.data$name),arid$name),X25.=rnorm(6,mean(arid$X25.,na.rm=TRUE),sd(arid$X25.,na.rm=TRUE)),X10.=rnorm(6,mean(arid$X10.,na.rm=TRUE),sd(arid$X10.,na.rm=TRUE)))
arid<-dplyr::select(arid,name,X25.,X10.)
arid<-rbind(arid,addos)
HF.data<-left_join(HF.data,arid)

###calculate the mean values for FLS and flowering tiem (flowers open for each species)
HF.means<-HF.data %>% group_by(name)%>% summarise(meanFLS.func=mean(funct.fls,na.rm=TRUE),meanFLS.phys=mean(phys.fls,na.rm=TRUE),meanFLS.inter=mean(inter.fls,na.rm=TRUE),meanflo=mean(flo_cent,na.rm=TRUE),meandrought=mean(precip_cent),meanpol=mean(pol_cent))
HF.means$name

####make them cateforical
HF.means$hyster.fuct<-ifelse(HF.means$meanFLS.func>1,1,0)
HF.means$hyster.phys<-ifelse(HF.means$meanFLS.phys>1,1,0)
HF.means$hyster.inter<-ifelse(HF.means$meanFLS.inter>1,1,0)

colnames(HF.data)
addmeans<-dplyr::select(HF.data,name,min_temp,seed_pound,moisture_use)

addmeans<-distinct(addmeans)
HF.means<-left_join(HF.means,addmeans,by="name")

HF.means$dummydrought=NA
HF.means$dummydrought[which(HF.means$moisture_use=="H")]<-1
HF.means$dummydrought[which(HF.means$moisture_use=="M")]<-.5
HF.means$dummydrought[which(HF.means$moisture_use=="L")]<-0
HF.means$mosit_cent<-(HF.means$dummydrought-mean(HF.means$dummydrought,na.rm=TRUE))/(2*sd(HF.means$dummydrought,na.rm=TRUE))
HF.means$mosit_cent<-ifelse(is.na(HF.means$mosit_cent),mean(HF.means$mosit_cent,na.rm=TRUE),HF.means$mosit_cent)

HF.means$cold_cent<-(HF.means$min_temp-mean(HF.means$min_temp,na.rm=TRUE))/(2*sd(HF.means$min_temp,na.rm=TRUE))
HF.means$seed_cent<-(HF.means$seed_pound-mean(HF.means$seed_pound,na.rm=TRUE))/(2*sd(HF.means$seed_pound,na.rm=TRUE))
HF.means$seed_cent<-ifelse(is.na(HF.means$seed_cent),mean(HF.means$seed_cent,na.rm=TRUE),HF.means$seed_cent)
HF.means$seed_cent<-(-(HF.means$seed_cent))

mtsv<-read.csv("MTSV_USFS/michdata_final.csv")  
fruiting<-data.frame(name=mtsv$name,fruiting=mtsv$fruiting)
HF.means<-left_join(HF.means,fruiting) ### add mean
HF.means$fruiting<-ifelse(is.na(HF.means$fruiting),mean(HF.means$fruiting,na.rm=TRUE),HF.means$fruiting)
HF.means$disperse_cent<-(HF.means$fruiting-mean(HF.means$fruiting,na.rm=TRUE))/(2*sd(HF.means$fruiting,na.rm=TRUE))

HF.means<-left_join(HF.means,arid)
HF.means$arid10_cent<-(HF.means$X10.-mean(HF.means$X10.,na.rm=TRUE))/(2*sd(HF.means$X10.,na.rm=TRUE))

######categorical models

HF.means<-  HF.means %>% remove_rownames %>% column_to_rownames(var="name")
z.phys.HFbin<-phylolm(hyster.phys~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE) ###fbb-lbb


z.phys.HFbin.cold<-phylolm(hyster.phys~meanflo+cold_cent+meanpol+meanflo:cold_cent+meanflo:meanpol+meanpol:cold_cent,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

z.phys.HFbin.moist<-phylolm(hyster.phys~meanflo+mosit_cent+meanpol+meanflo:mosit_cent+meanflo:meanpol+meanpol:mosit_cent,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                           start.beta=NULL, start.alpha=NULL,
                           boot=599,full.matrix = TRUE)


z.funct.HFbin<-phylolm(hyster.fuct~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE) ###fopn-l75


z.phys.HFbin.disperse<-phylolm(hyster.phys~disperse_cent+meandrought+meanpol+disperse_cent:meandrought+disperse_cent:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE) ###fbb-lbb

z.phys.HFbin.seedmass<-phylolm(hyster.phys~seed_cent+meandrought+meanpol+seed_cent:meandrought+seed_cent:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                               start.beta=NULL, start.alpha=NULL,
                               boot=599,full.matrix = TRUE) ###fbb-lbb


z.phys.HFbin.arid10<-phylolm(hyster.phys~meanflo+arid10_cent+meanpol+meanflo:arid10_cent+meanflo:meanpol+meanpol:arid10_cent,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE) ###fbb-lbb

##### quantitative means ()

z.funct.HFmean<-phylolm(meanFLS.func~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean) ###fbb-lbb

z.phys.HFmean<-phylolm(meanFLS.phys~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean) ###fopn-l75



###alternative models phys
z.phys.HFmean.cold<-phylolm(meanFLS.phys~meanflo+cold_cent+meanpol+meanflo:cold_cent+meanflo:meanpol+meanpol:cold_cent,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean) ###fopn-l75

z.phys.HFmean.disperse<-phylolm(meanFLS.phys~disperse_cent+meandrought+meanpol+disperse_cent:meandrought+disperse_cent:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.phys.HFmean.disperse)

z.phys.HFmean.moist<-phylolm(meanFLS.phys~meanflo+mosit_cent+meanpol+meanflo:mosit_cent+meanflo:meanpol+meanpol:mosit_cent,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)

z.phys.HFmean.seedmass<-phylolm(meanFLS.phys~seed_cent+meandrought+meanpol+seed_cent:meandrought+seed_cent:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.phys.HFmean.seedmass) ###fopn-l75

z.phys.HFmean.arid10<-phylolm(meanFLS.phys~meanflo+arid10_cent+meanpol+meanflo:arid10_cent+meanflo:meanpol+meanpol:arid10_cent,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
 ###fopn-l75

AIC(z.phys.HFmean)
AIC(z.phys.HFmean.cold)
AIC(z.phys.HFmean.moist)
AIC(z.phys.HFmean.disperse)
AIC(z.phys.HFmean.seedmass)
#z.inter.HFmean<-phylolm(meanFLS.inter~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
#summary(z.inter.HFmean)

#z.inter.HFbin<-phylolm(hyster.inter~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
 #                      start.beta=NULL, start.alpha=NULL,
  #boot=599,full.matrix = TRUE)


###Hierarchical models
inv.phylo <- MCMCglmm::inverseA(HF.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


modelcont.funct.wspecies.ind<-brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                  family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 

modelcont.phys.wspecies.ind<-brm(phys.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

#modelcont.inter.wspecies.ind<-brm(inter.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
 #                                family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)


##plotting  functions ########
summary(z.funct.HFbin)
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
}
extract_CIs<-function(x){
  dplyr::filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
}
extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}
#####################
#Plotting prep
#####################
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

####alternative drought plotting quant)
main.bin<-full_join(extract_coefs(z.phys.HFbin),extract_CIs(z.phys.HFbin),by="trait")
coldbin<-full_join(extract_coefs(z.phys.HFbin.cold),extract_CIs(z.phys.HFbin.cold),by="trait")
moistbin<-full_join(extract_coefs(z.phys.HFbin.moist),extract_CIs(z.phys.HFbin.moist),by="trait")
dispbin<-full_join(extract_coefs(z.phys.HFbin.disperse),extract_CIs(z.phys.HFbin.disperse),by="trait")
seedbin<-full_join(extract_coefs(z.phys.HFbin.seedmass),extract_CIs(z.phys.HFbin.seedmass),by="trait")
arid10<-full_join(extract_coefs(z.phys.HFbin.arid10),extract_CIs(z.phys.HFbin.arid10),by="trait")
AIC(z.phys.HFbin)
AIC(z.phys.HFbin.cold)
AIC(z.phys.HFbin.moist)

main.bin$model<-"main model"
coldbin$model<-"cold tolerance"
moistbin$model<-"moisture use"
dispbin$model<-"dispersal time"
seedbin$model<-"seed mass"
arid10$model<-"aridity index"
altdrought.bin<-rbind(main.bin,coldbin,moistbin,dispbin,seedbin,arid10)
colnames(altdrought.bin)<-c("trait","estimate","low","high","model")


prec<-full_join(extract_coefs(z.phys.HFmean),extract_CIs(z.phys.HFmean),by="trait")
col<-full_join(extract_coefs(z.phys.HFmean.cold),extract_CIs(z.phys.HFmean.cold),by="trait")
mois<-full_join(extract_coefs(z.phys.HFmean.moist),extract_CIs(z.phys.HFmean.moist),by="trait")
dispe<-full_join(extract_coefs(z.phys.HFmean.disperse),extract_CIs(z.phys.HFmean.disperse),by="trait")
seedy<-full_join(extract_coefs(z.phys.HFmean.seedmass),extract_CIs(z.phys.HFmean.seedmass),by="trait")
aridy<-full_join(extract_coefs(z.phys.HFmean.arid10),extract_CIs(z.phys.HFmean.arid10),by="trait")

prec$model<-"main model"
col$model<-"cold tolerance"
mois$model<-"moisture use"
dispe$model<-"dispersal time"
seedy$model<-"seed mass"
aridy$model<-"aridity index"
altdrought<-rbind(prec,col,mois,dispe,seedy,aridy)
colnames(altdrought)<-c("trait","estimate","low","high","model")

altdrought<-dplyr::filter(altdrought,trait!="(Intercept)")
altdrought<-dplyr::filter(altdrought,trait!="sigma2")
unique(altdrought$trait)
altdrought$trait[which(altdrought$trait=="meanpol")]<-"pollination syndrome"
altdrought$trait[which(altdrought$trait=="meanflo")]<- "early flowering"
altdrought$trait[which(altdrought$trait=="disperse_cent")]<- "early flowering"
altdrought$trait[which(altdrought$trait=="seed_cent")]<- "early flowering"

altdrought$trait[which(altdrought$trait=="meandrought")]  <- "water dynamics"
altdrought$trait[which(altdrought$trait=="cold_cent")]  <- "water dynamics"
altdrought$trait[which(altdrought$trait=="mosit_cent")]  <- "water dynamics"
altdrought$trait[which(altdrought$trait=="arid10_cent")]  <- "water dynamics"

altdrought$trait[which(altdrought$trait=="meandrought:meanpol")]<- "pollination:water dynamics"
altdrought$trait[which(altdrought$trait=="cold_cent:meanpol")]<- "pollination:water dynamics"
altdrought$trait[which(altdrought$trait=="mosit_cent:meanpol")]<- "pollination:water dynamics"
altdrought$trait[which(altdrought$trait=="arid10_cent:meanpol")]<- "pollination:water dynamics"

altdrought$trait[which(altdrought$trait=="meanflo:meanpol")]<-"pollination:early flowering"
altdrought$trait[which(altdrought$trait=="disperse_cent:meanpol")]<-"pollination:early flowering"
altdrought$trait[which(altdrought$trait=="seed_cent:meanpol")]<-"pollination:early flowering"


altdrought$trait[which(altdrought$trait=="meanflo:meandrought")]<-"early flowering:water dynamics"
altdrought$trait[which(altdrought$trait=="meanflo:cold_cent")]<-"early flowering:water dynamics"
altdrought$trait[which(altdrought$trait=="meanflo:mosit_cent")]<-"early flowering:water dynamics"
altdrought$trait[which(altdrought$trait=="meanflo:arid10_cent")]<-"early flowering:water dynamics"
altdrought$trait[which(altdrought$trait=="disperse_cent:meandrought")]<-"early flowering:water dynamics"
altdrought$trait[which(altdrought$trait=="seed_cent:meandrought")]<-"early flowering:water dynamics"
pd=position_dodgev(height=0.4)
alta<-altdrought %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(color=model),position=pd,size=2)+
  geom_errorbarh(aes(xmin=low,xmax=high,color=model),position=pd,height=0,size=.5)+
  ggthemes::theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flower and leaf budburst")+
  guides(size = "legend", linetype= "none")+theme(plot.margin = margin(r = 1.5, l = 1.5,t=20,b=10))+scale_color_brewer(type = "qual",palette = 7)+
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
     



altdrought.bin<-dplyr::filter(altdrought.bin,trait!="(Intercept)")
altdrought.bin<-dplyr::filter(altdrought.bin,trait!="sigma2")
unique(altdrought.bin$trait)
altdrought.bin$trait[which(altdrought.bin$trait=="meanpol")]<-"pollination syndrome"
altdrought.bin$trait[which(altdrought.bin$trait=="meanflo")]<- "early flowering"

altdrought.bin$trait[which(altdrought.bin$trait=="disperse_cent")]<- "early flowering"
altdrought.bin$trait[which(altdrought.bin$trait=="seed_cent")]<- "early flowering"

altdrought.bin$trait[which(altdrought.bin$trait=="meandrought")]  <- "water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="cold_cent")]  <- "water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="mosit_cent")]  <- "water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="arid10_cent")]  <- "water dynamics"

altdrought.bin$trait[which(altdrought.bin$trait=="meandrought:meanpol")]<- "pollination:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="cold_cent:meanpol")]<- "pollination:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="mosit_cent:meanpol")]<- "pollination:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="arid10_cent:meanpol")]<- "pollination:water dynamics"

altdrought.bin$trait[which(altdrought.bin$trait=="meanflo:meanpol")]<-"pollination:early flowering"
altdrought.bin$trait[which(altdrought.bin$trait=="disperse_cent:meanpol")]<-"pollination:early flowering"
altdrought.bin$trait[which(altdrought.bin$trait=="seed_cent:meanpol")]<-"pollination:early flowering"


altdrought.bin$trait[which(altdrought.bin$trait=="meanflo:meandrought")]<-"early flowering:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="meanflo:cold_cent")]<-"early flowering:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="meanflo:mosit_cent")]<-"early flowering:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="meanflo:arid10_cent")]<-"early flowering:water dynamics"

altdrought.bin$trait[which(altdrought.bin$trait=="disperse_cent:meandrought")]<-"early flowering:water dynamics"
altdrought.bin$trait[which(altdrought.bin$trait=="seed_cent:meandrought")]<-"early flowering:water dynamics"
pd=position_dodgev(height=0.4)
altb<-altdrought.bin %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(color=model),position=pd,size=2)+
  geom_errorbarh(aes(xmin=low,xmax=high,color=model),position=pd,height=0,size=.5)+
  ggthemes::theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flower and leaf budburst")+
  guides(size = "legend", linetype= "none")+theme(plot.margin = margin(r = 1.5, l = 1.5,t=20,b=10))+scale_color_brewer(type = "qual",palette = 7)


setEPS()
postscript("alternatepredictors.eps",width = 10, height = 4)
ggpubr::ggarrange(altb,alta,ncol=2,nrow=1,widths = c(2,1.2),labels=c("a)","b)"),hjust=c(-11.2,-.5),common.legend = TRUE, legend="right")
dev.off()
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
#########################################################

a<-binary.main %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(),position=pd,size=2,color="firebrick1")+
  geom_errorbarh(aes(xmin=low,xmax=high),position=pd,height=0,size=.5,color="firebrick1")+
  xlim(-1.5,2)+ ggthemes::theme_base(base_size = 10)+geom_vline(aes(xintercept=0),color="black")+xlab("Likelihood of hysteranthy")+
 guides(size = "legend", linetype= "none")+theme(plot.margin = margin(l=5,r = 1,t=20,b=10))+
  annotate("text",x=-1.3,y=6, label="biotic",fontface =3, size=2.5)+
  annotate("text",x=1.8,y=6, label="wind",fontface =3, size=2.5)+
  annotate("text",x=-1.3,y=5, label="drier",fontface =3, size=2.5)+
  annotate("text",x=1.8,y=5, label="wetter",fontface =3, size=2.5)+
  annotate("text",x=-1.3,y=4, label="earlier",fontface =3, size=2.5)+
  annotate("text",x=1.8,y=4, label="later",fontface =3, size=2.5) #This is the binary plots for fbb-lbb
  

b<-quant.main %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(),position=pd,size=2,color="goldenrod3")+
  geom_errorbarh(aes(xmin=low,xmax=high),position=pd,height=0,size=.5,color="goldenrod3")+
  xlim(-40,50)+ggthemes::theme_base(base_size = 10)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flower and leaf budburst")+
 guides(size = "legend", linetype= "none")+theme(axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                                                plot.margin = margin(r = 1.5, l = 1.5,t=20,b=10))+
  annotate("text",x=-36,y=6, label="biotic",fontface =3, size=2.5)+
  annotate("text",x=47,y=6, label="wind",fontface =3, size=2.5)+
  annotate("text",x=-36,y=5, label="drier",fontface =3, size=2.5)+
  annotate("text",x=47,y=5, label="wetter",fontface =3, size=2.5)+
  annotate("text",x=-36,y=4, label="earlier",fontface =3, size=2.5)+
  annotate("text",x=47,y=4, label="later",fontface =3, size=2.5)
                                                                

cont.phys<-filter(cont,class=="fbb-lbb")
c<-cont.phys %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(),position=pd,size=2,stroke=.5,color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,height=0,linetype="dotted",color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,height=0,linetype="solid",color="darkorchid3")+
  ggthemes::theme_base(base_size = 10)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flower and leaf budburst")+
  xlim(-40,50)+theme(axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                                    plot.margin = margin(l = 1,r=5,t=20,b=10) )+
  annotate("text",x=-36,y=6, label="biotic",fontface =3, size=2.5)+
  annotate("text",x=46,y=6, label="wind",fontface =3, size=2.5)+
  annotate("text",x=-36,y=5, label="drier",fontface =3, size=2.5)+
  annotate("text",x=46,y=5, label="wetter",fontface =3, size=2.5)+
  annotate("text",x=-36,y=4, label="earlier",fontface =3, size=2.5)+
  annotate("text",x=46,y=4, label="later",fontface =3, size=2.5)
  



setEPS()
postscript("HFmodelplots.eps",width = 10, height = 4)
ggpubr::ggarrange(a,b,c,nrow=1,ncol=3,widths=c(1.7,1,1),labels=c("a)","b)","c)"),hjust=c(-11.2,-.5,-.5))
dev.off()

#tiff("HFmodplots.tiff",width = 8, height = 4,units = "in",res= 200)
#ggpubr::ggarrange(a,b,c,nrow=1,ncol=3,widths=c(2,1,1),labels=c("a)","b)","c)"),hjust=c(-11.2,-.5,-.5))
#dev.off()
######Suppliment

bin.comps<-filter(comps,model=="categorical")
quant.comps<-filter(comps,model!="categorical")

supa<-bin.comps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class),position=pd,size=2,color="firebrick1")+
  geom_errorbarh(aes(xmin=low,xmax=high,group=class),position=pd,height=0,size=.5,color="firebrick1")+
  xlim(-1.5,2)+ ggthemes::theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlab("Likelihood of hysteranthy")+
  guides(size = "legend", linetype= "none")+theme(
                                                  plot.margin = margin(l=5,r = 1,t=20,b=10))+
  annotate("text",x=-1.3,y=6, label="biotic",fontface =3, size=2.5)+
  annotate("text",x=1.8,y=6, label="wind",fontface =3, size=2.5)+
  annotate("text",x=-1.3,y=5, label="drier",fontface =3, size=2.5)+
  annotate("text",x=1.8,y=5, label="wetter",fontface =3, size=2.5)+
  annotate("text",x=-1.3,y=4, label="earlier",fontface =3, size=2.5)+
  annotate("text",x=1.8,y=4, label="later",fontface =3, size=2.5) #This is the binary plots for fbb-lbb


supb<-quant.comps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class),position=pd,size=2,color="goldenrod3")+
  geom_errorbarh(aes(xmin=low,xmax=high,group=class),position=pd,height=0,size=.5,color="goldenrod3")+
  xlim(-45,50)+ggthemes::theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flowering and leafing")+
  guides(size = "legend", linetype= "none")+theme(axis.title.y=element_blank(),
                                                  axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                                  plot.margin = margin(r = 1.5, l = 1.5,t=20,b=10))+
  
  annotate("text",x=-36,y=6, label="biotic",fontface =3, size=2.5)+
  annotate("text",x=46,y=6, label="wind",fontface =3, size=2.5)+
  annotate("text",x=-36,y=5, label="drier",fontface =3, size=2.5)+
  annotate("text",x=46,y=5, label="wetter",fontface =3, size=2.5)+
  annotate("text",x=-36,y=4, label="earlier",fontface =3, size=2.5)+
  annotate("text",x=46,y=4, label="later",fontface =3, size=2.5)

supc<-cont %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=class),position=pd,size=2,stroke=.5,color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=class),position=pd,height=0,linetype="dotted",color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,group=class),position=pd,height=0,linetype="solid",color="darkorchid3")+
  ggthemes::theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flowering and leafing")+
  xlim(-45,50)+theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                    
                     axis.ticks.y=element_blank(),
                     plot.margin = margin(l = 1,r=5,t=20,b=10) )+
   annotate("text",x=-36,y=6, label="biotic",fontface =3, size=2.5)+
  annotate("text",x=46,y=6, label="wind",fontface =3, size=2.5)+
  annotate("text",x=-36,y=5, label="drier",fontface =3, size=2.5)+
  annotate("text",x=46,y=5, label="wetter",fontface =3, size=2.5)+
  annotate("text",x=-36,y=4, label="earlier",fontface =3, size=2.5)+
  annotate("text",x=46,y=4, label="later",fontface =3, size=2.5)

setEPS()
postscript("HFmodelplots4SUPP.eps",width = 8, height = 4)
ggpubr::ggarrange(supa,supb,supc,nrow=1,ncol=3,widths=c(2,1,1),labels=c("a)","b)","c)"),hjust=c(-11.2,-.5,-.5),common.legend = TRUE,legend = "right")

dev.off()
data.frame(phenophase=c("leaf budburst","flower budburst","flowers open", "75% of leaves full size"),HF=c("lbb","fbb","fopn","l75"),BBCH=c("09","55","60","17"))

##### alternative models
colnames(HF.data)


modelcont.phys.wspecies.ind<-brm(phys.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

###### alternative models:



HF.data$dummydrought=NA
HF.data$dummydrought[which(HF.data$moisture_use=="H")]<-1
HF.data$dummydrought[which(HF.data$moisture_use=="M")]<-.5
HF.data$dummydrought[which(HF.data$moisture_use=="L")]<-0
HF.data$dummy_cent<-(HF.data$dummydrought-mean(HF.data$dummydrought,na.rm=TRUE))/(2*sd(HF.data$dummydrought,na.rm=TRUE))


modelcont.phys.moist<-brm(phys.fls~ pol+flo_cent+dummy_cent+dummy_cent:flo_cent+dummy_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

moist.cont<-extract_coefs4HF(modelcont.phys.moist)

xtable(phys.cont,digits=3,label = "HF_phys_main",caption = "Model results for main text hierarchical model")
xtable(moist.cont,digits=3,label = "HF_phys_moisture",caption = "Model results for alternative hierarchical model with moisture use as a proxy for the water limitation hypothesis")

  ggplot(moist.cont,aes(Estimate,trait))+geom_point(aes(),position=pd,size=2,stroke=.5,color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5),position=pd,height=0,linetype="dotted",color="darkorchid3")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90),position=pd,height=0,linetype="solid",color="darkorchid3")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  xlim(-35,50)

  cor(HF.data$dummydrought,HF.data$min_precip,use = "pairwise.complete.obs")
?cor()
  

HF.data$disperse_cent<-(HF.data$fruiting-mean(HF.data$fruiting,na.rm=TRUE))/(2*sd(HF.data$fruiting,na.rm=TRUE))
modelcont.phys.disperse.moist<-brm(phys.fls~ pol+disperse_cent+precip_cent+precip_cent:disperse_cent+precip_cent:pol+pol:disperse_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99),iter=8000, warmup=6000)

###seed mass HF
addins<-read.csv("mich_tree_additions.csv")
HF.data<-dplyr::left_join(HF.data,addins,by="name")
table(HF.data$seed_mass)
table(HF.data$seed_pound)
HF.data$seedmass_centy<-(HF.data$seed_pound-mean(HF.data$seed_pound,na.rm=TRUE))/(2*sd(HF.data$seed_pound,na.rm=TRUE))
HF.data$seedmass_cent<--(HF.data$seedmass_centy)## take the inverse

modelcont.phys.seedmass<-brm(phys.fls~ pol+seedmass_cent+precip_cent+precip_cent:seedmass_cent+precip_cent:pol+pol:seedmass_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                   family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99,max_treedepth=12),iter=8000, warmup=7000)




#summary(modelcont.phys.disperse.moist)

#moist.disp.cont<-extract_coefs4HF(modelcont.phys.disperse.moist)

#xtable(moist.disp.cont,digits=3,label = "HF_phys_moist_disp",caption = "Model results for alternative hierarchical model with moisture use as a proxy for the water limitation and dispersal time as a proxy for the early flowering hypothesis")



HF.data$arid.25.cent<-(HF.data$X25.-mean(HF.data$X25.,na.rm=TRUE))/(2*sd(HF.data$X25.,na.rm=TRUE))
HF.data$arid.10.cent<-(HF.data$X10.-mean(HF.data$X10.,na.rm=TRUE))/(2*sd(HF.data$X10.,na.rm=TRUE))

modelcont.phys.wspecies.aridity<-brm(phys.fls~ pol+flo_cent+arid.10.cent+arid.10.cent:flo_cent+arid.10.cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                     family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.99,max_treedepth=20),iter=6500, warmup=5000)



save.image("finalmanuscriptmodels") 

###checkk PETP
climate<-read.csv("PETP.HF.csv")
HF.data<-left_join(HF.data,climate)
HF.data$aridindex_cent<-(HF.data$aridindex-mean(HF.data$aridindex,na.rm=TRUE))/(2*sd(HF.data$aridindex,na.rm=TRUE))
HF.data$PETP_cent<-(HF.data$PTEP-mean(HF.data$PTEP,na.rm=TRUE))/(2*sd(HF.data$PTEP,na.rm=TRUE))

PETP.mod<-brm(phys.fls~PETP_cent*flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                  family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

summary(PETP.mod)
arid.mod<-brm(phys.fls~aridindex_cent*flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
              family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

full.petep<-brm(phys.fls~PETP_cent+flo_cent+PETP_cent:pol_cent+pol_cent:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

summary(full.petep)
cor(HF.data$flo_cent,HF.data$PETP_cent,use = "pairwise.complete.obs" )



just.PETP.mod<-brm(phys.fls~PTEP+(1|name)+(1|tree.id/species), data = HF.data, 
              family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)
extract_coefs4HF(just.PETP.mod)

syn.PETP.mod<-brm(phys.fls~PTEP+PTEP:pol+(1|name)+(1|tree.id/species), data = HF.data, 
                   family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)


extract_coefs4HF(full.petep)

summary(just.PETP.mod)
just.arid.mod<-brm(phys.fls~aridindex+(1|name)+(1|tree.id/species), data = HF.data, 
              family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000) 


PETP<-extract_coefs4HF(just.PETP.mod)
xtable(PETP,digits=3,label = "HF.PETP",caption = "Model results from Bayesian hierarchical regresssion analysis 
       show that annual variation in water balance (P-ETP) does not significantly influence FLS variation at Harvard Forest")



### alternative models with cold tolerance

cor(HF.data$min_precip,HF.data$min_temp,use = "pairwise.complete.obs")
cor(HF.data$dummydrought,HF.data$min_temp,use = "pairwise.complete.obs")
cor(HF.data$dummydrought,HF.data$min_precip,use = "pairwise.complete.obs")

cor(HF.means$X50.,HF.means$meandrought,use = "pairwise.complete.obs")

HF.data$cold_cent<-(HF.data$min_temp-mean(HF.data$min_temp))/(2*sd(HF.data$min_temp))

modelcont.phys.wcols.ind<-brm(phys.fls~ pol+flo_cent+cold_cent+cold_cent:flo_cent+cold_cent:pol+pol:flo_cent+(1|name)+(1|tree.id/species), data = HF.data, 
                                 family = gaussian(), cov_ranef = list(name= A),control=list(adapt_delta=0.95),iter=4000, warmup=3000)

coldy<-extract_coefs4HF(modelcont.phys.wcols.ind)
xtable(coldy,digits=3,label = "HF.coldtol",caption = "Cold tolerance instead of water limitation")

###plots for alternative models
coldy$model<-"min. T"
moist.cont$model<-"moisture use"
cont.phys$model<-"min. P"
cont.phys<-dplyr::select(cont.phys,-class)

alternateuniverse<-rbind(coldy,moist.cont,cont.phys)
alternateuniverse<-dplyr::filter(alternateuniverse,trait!="(Intercept)")
alternateuniverse<-dplyr::filter(alternateuniverse,trait!="Intercept")
#comps<-dplyr::filter(comps,trait!="sigma2")
unique(alternateuniverse$trait)
alternateuniverse$trait[which(alternateuniverse$trait=="pol")]<-"pollination syndrome"
alternateuniverse$trait[which(alternateuniverse$trait=="flo_cent")]<- "early flowering"
alternateuniverse$trait[which(alternateuniverse$trait=="cold_cent")]  <- "water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="dummy_cent")]  <- "water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="water dynamics")]  <- "water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="pol:flo_cent")]<-"pollination:early flowering"

alternateuniverse$trait[which(alternateuniverse$trait=="flo_cent:dummy_cent")]<-"early flowering:water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="flo_cent:cold_cent")]<-"early flowering:water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="earlier flowering:water dynamics")]<-"early flowering:water limitation"

alternateuniverse$trait[which(alternateuniverse$trait=="earlier flowering")]<- "early flowering"

alternateuniverse$trait[which(alternateuniverse$trait=="pol:dummy_cent")]<-"pollination:water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="pol:cold_cent")]<-"pollination:water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="pollination:water dynamics")]<-"pollination:water limitation"
alternateuniverse$trait[which(alternateuniverse$trait=="pollination:earlier flowering")]<-"pollination:early flowering"

pd=position_dodgev(height=0.4)

#alta<alternateuniverse %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water limitation","pollination:early flowering","pollination:water limitation","early flowering","water limitation","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(color=model),position=pd,size=2,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=model),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=model),position=pd,height=0,linetype="solid")+
  theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+scale_color_brewer(type="qual" ,palette = 2 )+
  xlim(-45,50)+xlab("Days between flower and leaf budburst")

#####
disp.cont<-extract_coefs4HF(modelcont.phys.disperse.moist)
disp.cont$model<-"dispersal time"

seed.cont<-extract_coefs4HF(modelcont.phys.seedmass)
seed.cont$model<-"seeds mass"
contphys2<-dplyr::select(cont.phys,-class)
contphys2$model= "flowering time"
alternateflo<-rbind(disp.cont,seed.cont,contphys2)
alternateflo<-dplyr::filter(alternateflo,trait!="Intercept")
#comps<-dplyr::filter(comps,trait!="sigma2")
unique(alternateflo$trait)
alternateflo$trait[which(alternateflo$trait=="pol")]<-"pollination syndrome"
alternateflo$trait[which(alternateflo$trait=="flo_cent")]<- "early flowering"
alternateflo$trait[which(alternateflo$trait=="disperse_cent")]<- "early flowering"
alternateflo$trait[which(alternateflo$trait=="seedmass_cent")]<- "early flowering"
alternateflo$trait[which(alternateflo$trait=="earlier flowering")]<- "early flowering"
alternateflo$trait[which(alternateflo$trait=="water dynamics")]  <- "water limitation"
alternateflo$trait[which(alternateflo$trait=="precip_cent")]  <- "water limitation"
alternateflo$trait[which(alternateflo$trait=="pol:disperse_cent")]<-"pollination:early flowering"

alternateflo$trait[which(alternateflo$trait=="earlier flowering:water dynamics")]<-"early flowering:water limitation"
alternateflo$trait[which(alternateflo$trait=="pollination:water dynamics")]<-"pollination:water limitation"
alternateflo$trait[which(alternateflo$trait=="pollination:earlier flowering")]<-"pollination:early flowering"

alternateflo$trait[which(alternateflo$trait=="pol:precip_cent")]<-"pollination:water limitation"
alternateflo$trait[which(alternateflo$trait=="disperse_cent:precip_cent")]<-"early flowering:water limitation"
alternateflo$trait[which(alternateflo$trait=="seedmass_cent:precip_cent")]<-"early flowering:water limitation"
alternateflo$trait[which(alternateflo$trait=="pol:disperse_cent")]<-"pollination:early flowering"
alternateflo$trait[which(alternateflo$trait=="pol:seedmass_cent")]<-"pollination:early flowering"
alternateflo$model[which(alternateflo$model=="min. P")]<-"flowering time"

altb<-alternateflo %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water limitation","pollination:early flowering","pollination:water limitation","early flowering","water limitation","pollination syndrome"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(color=model),position=pd,size=2,stroke=.5)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=model),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=model),position=pd,height=0,linetype="solid")+scale_color_brewer(type="qual" ,palette = 3 )+
  ggthemes::theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlab("Days between flower and leaf budburst")


ggpubr::ggarrange(alta,altb)

(-2.806514e-16)*(2*sd(HF.data$fbb.jd,na.rm=TRUE))+mean(HF.data$fbb.jd,na.rm=TRUE)

#save.image("altpredictors") 
######3marginal effects
mean(HF.data$flo_cent,na.rm=TRUE)
apc.funct<- ggeffects::ggpredict(modelcont.phys.wspecies.ind,c("precip_cent","pol","flo_cent[-2.806514e-16]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot<-plot(apc.funct)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1,1.5),labels=c(6,13,19,26,33,40,47))+
  xlab("Min. precipitation across range (cm)")+ylab("Flower to leaf budburst (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag=NULL)+ggthemes::theme_base(base_size = 11)

mean(HF.data$min_precip)
range(HF.data$precip_cent)
options(scipen = 999)
((mean(HF.data$min_precip))*(2*sd(HF.data$precip_cent,na.rm=TRUE)))+mean(HF.data$precip_cent,na.rm=TRUE)


(.5)*(2*sd(HF.data$fbb.jd,na.rm=TRUE))+mean(HF.data$fbb.jd,na.rm=TRUE)

apc.funct2<- ggeffects::ggpredict(modelcont.phys.wspecies.ind,c("flo_cent","pol","precip_cent[-2.253749e-17]"), ci.lvl=0.50)  #May the fourth
apc.funct.plot2<-plot(apc.funct2)+scale_x_continuous(breaks =c(-1.0,0,1),labels=c(84,128,171))+
  xlab("Flowering time (day of year)")+ylab("Flower to leaf budburst (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag=NULL)+ggthemes::theme_base(base_size = 11)

apc.funct3<- ggeffects::ggpredict(modelcont.phys.wspecies.ind,c("flo_cent","precip_cent","pol[0]"), ci.lvl=0.50)
apc.funct.plot3<-plot(apc.funct3)+scale_x_continuous(breaks =c(-1.0,0,1),labels=c(84,128,171))+
  xlab("Flowering time (day of year)")+ylab("Flower to leaf budburst (days)")+scale_colour_manual(name="Minimum precipitation",labels=c("low","average","high"),values=c("red","yellow","blue"))+scale_fill_manual(name="minimum precipitation",labels=c("low","average","high"),values=c("red","yellow","blue"))+
  labs(title = NULL,tag=NULL)+ggthemes::theme_base(base_size = 11)


#tiff("apcs.tif",width = 8, height = 4,units = "in",res = 200 )
jpeg("apcs.jpeg",width = 8, height = 4,units = "in",res = 200 )
ggpubr::ggarrange(apc.funct.plot,apc.funct.plot2,nrow=1,common.legend = TRUE,legend="right",labels = c("a)","b)"))
dev.off()


apc.funct.plot
dev.off()


