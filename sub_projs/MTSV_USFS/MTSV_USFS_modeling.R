##modeling USFS and MTSV


rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_projs/MTSV_USFS/")

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
#load("MTSVUSFS.mods")

mich.data<-read.csv("michdata_final.csv")
mich.tre<-read.tree("michtre_final.tre")

silv.data<-read.csv("silvdata_final.csv")
silv.tre<-read.tree("silvtre_final.tre")
##### phylo signal for suppliment
mich.tre$node.label<-NULL
silv.tre$node.label<-NULL

##center predictors for 
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$precip_cent<-(mich.data$min._precip-mean(mich.data$min._precip))/(2*sd(mich.data$min._precip)) 

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
z.funct.MTSV<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=599,full.matrix = TRUE)


z.phys.MTSV<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

#z.flex.MTSV<-phyloglm(pro~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#                      start.beta=NULL, start.alpha=NULL,
 #                     boot=599,full.matrix = TRUE)

#USFS models
z.funct.USFS<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)


z.phys.USFS<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

#z.flex.USFS<-phyloglm(pro~pol_cent+flo_cent.neg+precip.neg+precip.neg:flo_cent.neg+precip.neg:pol_cent+pol_cent:flo_cent.neg,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
 #                     start.beta=NULL, start.alpha=NULL,
  #                    boot=599,full.matrix = TRUE)

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

#MTSV.flex.dat<-full_join(extract_coefs(z.flex.MTSV),extract_CIs(z.flex.MTSV),by="trait")
# colnames(MTSV.flex.dat)<-c("trait","estimate","low","high")
# MTSV.flex.dat$class<-"intermediate"

MTSV<-rbind(MTSV.funct.dat,MTSV.phys.dat)
MTSV$data<-"MTSV"   

##USFS
USFS.funct.dat<-full_join(extract_coefs(z.funct.USFS),extract_CIs(z.funct.USFS),by="trait")

colnames(USFS.funct.dat)<-c("trait","estimate","low","high")
USFS.funct.dat$class<-"functional"

USFS.phys.dat<-full_join(extract_coefs(z.phys.USFS),extract_CIs(z.phys.USFS),by="trait")
colnames(USFS.phys.dat)<-c("trait","estimate","low","high")
USFS.phys.dat$class<-"physiological"

#USFS.flex.dat<-full_join(extract_coefs(z.flex.USFS),extract_CIs(z.flex.USFS),by="trait")
#colnames(USFS.flex.dat)<-c("trait","estimate","low","high")
#USFS.flex.dat$class<-"intermediate"

USFS<-rbind(USFS.funct.dat,USFS.phys.dat)
USFS$data<-"USFS"

comps<-rbind(USFS,MTSV)
comps<-dplyr::filter(comps,trait!="(Intercept)")

comps$trait[which(comps$trait=="pol_cent")]<-"pollination syndrome"
comps$trait[which(comps$trait=="flo_cent")]<- "earlier flowering"
comps$trait[which(comps$trait=="precip_cent")]  <- "water dynamics"
comps$trait[which(comps$trait=="pol_cent:precip_cent")]<- "pollination:water dynamics"
comps$trait[which(comps$trait=="pol_cent:flo_cent")]<-"pollination:earlier flowering"
comps$trait[which(comps$trait=="flo_cent:precip_cent")]<-"earlier flowering:water dynamics"


jpeg("MTSV.USFS.jpeg",width = 8, height = 4, units = 'in', res=200)
pd=position_dodgev(height=0.4)
comps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("earlier flowering:water dynamics","pollination:earlier flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class,color=data),position=pd,size=3)+
  geom_errorbarh(aes(xmin=low,xmax=high,linetype=class,color=data),position=pd,height=0,size=1)+
 scale_linetype_manual(values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-8,9)+
  scale_color_manual(values=c("orchid4","springgreen4"))+guides(size = "legend", linetype= "none")
  
dev.off()

  save.image("MTSVUSFS.mods")  
  