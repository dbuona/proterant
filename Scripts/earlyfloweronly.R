####c### This is the final analysis file for hysteranthy anaylsis on MTSV as of 3/28/18.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
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

mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time


###restrict data set to early flowering
median(mich.data$flo_time)
mean(mich.data$flo_time)
earlymich<-filter(mich.data,flo_time<=5.5)

###prune tree
mich.data<-rownames_to_column(mich.data, "name")
namelist<-earlymich$name
names.intree<-mich.tree$tip.label

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
earlytree<-drop.tip(mich.tree,to.prune)

earlytree$tip.label==earlymich$name

###For some reason I need to runthis twice when sourcing
earlymich<-filter(mich.data,flo_time<=5.5)

###prune tree
#mich.data<-rownames_to_column(mich.data, "name")
namelist<-earlymich$name
names.intree<-mich.tree$tip.label

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
earlytree<-drop.tip(mich.tree,to.prune)

earlytree$tip.label==earlymich$name

#recenter
earlymich$height_cent<-(earlymich$heigh_height-mean(earlymich$heigh_height))/(2*sd(earlymich$heigh_height))
earlymich$fruit_cent<-(earlymich$fruiting-mean(earlymich$fruiting))/(2*sd(earlymich$fruiting))
earlymich$flo_cent<-(earlymich$flo_time-mean(earlymich$flo_time))/(2*sd(earlymich$flo_time))
earlymich$pol_cent<-(earlymich$pol-mean(earlymich$pol))/(2*sd(earlymich$pol))
earlymich$dev_time_center<-(earlymich$dev.time-mean(earlymich$dev.time))/(2*sd(earlymich$dev.time))

earlymich<-  earlymich %>% remove_rownames %>% column_to_rownames(var="name")

early.seed.cent<-phyloglm(pro2~pol+height_cent+flo_cent+dev_time_center+shade_bin,earlymich, earlytree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=599,full.matrix = TRUE)
summary(early.seed.cent)

bootest<-as.data.frame(early.seed.cent$coefficients)
bootconf<-as.data.frame(early.seed.cent$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="alpha")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
bootmich$trait[bootmich$trait=="shade_bin"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="dev_time_center"]<-"seed development"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

jpeg("early.effectplot.jpeg")
functplot1<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
functplot1
dev.off()

