####c### This is the final analysis file for hysteranthy anaylsis on MTSV as of 3/28/18.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
#graphics.off()
setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)
library("randomForest")

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
drought.dat<-read.csv("..//Data/USDA_traitfor_MTSV.csv",header=TRUE)

mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time
###one more cleaninging tax
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)

###make the tree work 
mich.tree$node.label<-NULL


######### later we'll be running model with drought tolerance, so add drought tolerance data
mich.data<-left_join(mich.data,drought.dat,by="name")

######prune the tree for drought modeling
mich.data.wdrought<-filter(mich.data,!is.na(min._precip))
mich.data.wdrought$precip_cent<-(mich.data.wdrought$min._precip-mean(mich.data.wdrought$min._precip))/(2*sd(mich.data.wdrought$min._precip))
earlymich<-filter(mich.data.wdrought,flo_time<=5.5)

names.intree<-mich.tree$tip.label
namelist<-unique(earlymich$name)
to.prune<-which(!names.intree%in%namelist)
mich.tree.early<-drop.tip(mich.tree,to.prune)
mytree.names<-mich.tree.early$tip.label

###Rescale predictors this makes it so you can compare binary to continous data


earlymich$flo_cent<-(earlymich$flo_time-mean(earlymich$flo_time))/(2*sd(earlymich$flo_time))
earlymich$pol_cent<-(earlymich$pol-mean(earlymich$pol))/(2*sd(earlymich$pol))
earlymich$precip_cent<-(earlymich$min._precip-mean(earlymich$min._precip))/(2*sd(earlymich$min._precip))




earlymich<-  earlymich %>% remove_rownames %>% column_to_rownames(var="name")

z.funct.early<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,earlymich, mich.tree.early, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

z.phys.early<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,earlymich, mich.tree.early, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)



bootestSF<-as.data.frame(z.funct.early$coefficients)
bootconfSF<-as.data.frame(z.funct.early$bootconfint95)
bootconfSF<-as.data.frame(t(bootconfSF))
bootestSF<-rownames_to_column(bootestSF, "trait")
bootconfSF<-rownames_to_column(bootconfSF, "trait")
bootdroughtF<-full_join(bootconfSF,bootestSF, by="trait")
colnames(bootdroughtF)<-c("trait","low","high","estimate")
bootdroughtF<-dplyr::filter(bootdroughtF, trait!="alpha")
bootdroughtF<-dplyr::filter(bootdroughtF, trait!="(Intercept)")
bootdroughtF$class<-"functional-MTSV"

bootestSP<-as.data.frame(z.phys.early$coefficients)
bootconfSP<-as.data.frame(z.phys.early$bootconfint95)
bootconfSP<-as.data.frame(t(bootconfSP))
bootestSP<-rownames_to_column(bootestSP, "trait")
bootconfSP<-rownames_to_column(bootconfSP, "trait")
bootdroughtP<-full_join(bootconfSP,bootestSP, by="trait")
colnames(bootdroughtP)<-c("trait","low","high","estimate")
bootdroughtP<-dplyr::filter(bootdroughtP, trait!="alpha")
bootdroughtP<-dplyr::filter(bootdroughtP, trait!="(Intercept)")
bootdroughtP$class<-"physiological-MTSV"

###combine these three catagories
bootearly<-rbind(bootdroughtF,bootdroughtP)
#bootdrought<-rbind(bootdrought,bootdroughtI)

#and plot
pd=position_dodgev(height=0.3)
plotty<-ggplot(bootdrought,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()
plotty



