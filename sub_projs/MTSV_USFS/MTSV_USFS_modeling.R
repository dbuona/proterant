##modeling USFS and MTSV


rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")

###libraryies
library("ggplot2")
library(brms)
library(ggthemes)
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library(boot)
library("phylolm")
library(ggstance)
library(tibble)
library(dplyr)

library("raster")
library("remote")
library(reshape2)
library(RColorBrewer)

mich.data<-read.csv("datasheets_derived/MTSV_USFS/michdata_final.csv")
mich.tre<-read.tree("datasheets_derived/MTSV_USFS/michtre_final.tre")

silv.data<-read.csv("datasheets_derived/MTSV_USFS/silvdata_final.csv")
silv.tre<-read.tree("datasheets_derived/MTSV_USFS/silvtre_final.tre")
##### phylo signal for suppliment
mich.tre$node.label<-NULL
silv.tre$node.label<-NULL

##center predictors for USFS

silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))
silv.data$pol_cent<-(silv.data$pol-mean(silv.data$pol))/(2*sd(silv.data$pol))
silv.data$precip_cent<-(silv.data$min._precip-mean(silv.data$min._precip))/(2*sd(silv.data$min._precip)) 

######### phyloglm requires species names to be in rownames
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<- silv.data %>% remove_rownames %>% column_to_rownames(var="name")

##flip axis for confusing negative predictors
mich.data$flo_cent.neg<--(mich.data$flo_cent)
silv.data$flo_cent.neg<--(silv.data$flo_cent)

mich.data$precip.neg<--(mich.data$precip_cent)
silv.data$precip.neg<--(silv.data$precip_cent)

###mtsv models
z.funct.MTSV<-phyloglm(pro2~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)


z.phys.MTSV<-phyloglm(pro3~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

z.flex.MTSV<-phyloglm(pro~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

#USFS models
z.funct.USFS<-phyloglm(pro2~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)


z.phys.USFS<-phyloglm(pro3~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

z.flex.USFS<-phyloglm(pro~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

###function for plotting
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
}
extract_CIs<-function(x){
  dplyr::filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
}

####pull MTSV together
MTSV.funct.dat<-full_join(extract_coefs(z.funct.MTSV),extract_CIs(z.funct.MTSV),by="trait")

colnames(MTSV.funct.dat)<-c("trait","estimate","low","high")
MTSV.funct.dat$class<-"functional"

MTSV.phys.dat<-full_join(extract_coefs(z.phys.MTSV),extract_CIs(z.phys.MTSV),by="trait")
colnames(MTSV.phys.dat)<-c("trait","estimate","low","high")
MTSV.phys.dat$class<-"physiological"

MTSV.flex.dat<-full_join(extract_coefs(z.flex.MTSV),extract_CIs(z.flex.MTSV),by="trait")
 colnames(MTSV.flex.dat)<-c("trait","estimate","low","high")
 MTSV.flex.dat$class<-"intermediate"

MTSV<-rbind(MTSV.funct.dat,MTSV.flex.dat,MTSV.phys.dat)
MTSV$data<-"MTSV"   

##USFS
USFS.funct.dat<-full_join(extract_coefs(z.funct.USFS),extract_CIs(z.funct.USFS),by="trait")

colnames(USFS.funct.dat)<-c("trait","estimate","low","high")
USFS.funct.dat$class<-"functional"

USFS.phys.dat<-full_join(extract_coefs(z.phys.USFS),extract_CIs(z.phys.USFS),by="trait")
colnames(USFS.phys.dat)<-c("trait","estimate","low","high")
USFS.phys.dat$class<-"physiological"

USFS.flex.dat<-full_join(extract_coefs(z.flex.USFS),extract_CIs(z.flex.USFS),by="trait")
colnames(USFS.flex.dat)<-c("trait","estimate","low","high")
USFS.flex.dat$class<-"intermediate"

USFS<-rbind(USFS.funct.dat,USFS.flex.dat,MTSV.phys.dat)
USFS$data<-"USFS"

comps<-dplyr::filter(comps,trait!="(Intercept)")

comps$trait[which(comps$trait=="pol_cent")]<-"pollination syndrome"
comps$trait[which(comps$trait=="flo_cent.neg")]<- "earlier flowering"
comps$trait[which(comps$trait=="precip.neg")]  <- "water dynamics"
comps$trait[which(comps$trait=="pol_cent:precip.neg")]<- "pollination:water dynamics"
comps$trait[which(comps$trait=="pol_cent:flo_cent.neg")]<-"pollination:flowering"
comps$trait[which(comps$trait=="flo_cent.neg:precip.neg")]<-"flowering:water dynamics"


pd=position_dodgev(height=0.4)
ggplot(comps,aes(estimate,trait))+geom_point(size=4,aes(color=class,shape=data),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class,linetype=data))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)


jpeg("..//sub_projs//MTSV.USFS.jpeg",width = 8.6, height = 4, units = 'in', res=200)
pd=position_dodgev(height=0.4)
comps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=data,color=class),position=pd,size=3)+
  geom_errorbarh(aes(xmin=low,xmax=high,linetype=data,color=class),position=pd,height=0)+
 scale_linetype_manual(values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-10,10)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+guides(size = "legend", linetype= "none")+
  annotate("text", x = 8.2, y = 6.4, label = "Hysteranthy",fontface="bold",size=3)+annotate("text", x = -8.9, y = 6.4, label = "Seranthy",fontface="bold",size=3)
  dev.off()
  
  