###This is a tidy analysis compiling and comparing all of Dan B's hysteranthy codes from elsewhere
###Began 5 Spet 2017

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(tidyverse)
library(boot)
library(phylolm)

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")

silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")

keeler.tree<-read.tree("pruned_keeler.tre")
keeler.data<-read.csv("keeler_cleaned.csv")

setdiff(keeler.data$name,mich.data$name)

######SET UP COMPARISON DATA SETS
names.intree<-mich.tree$tip.label
mich.by.keel<-intersect(keeler.data$name,mich.data$name)
to.prune<-which(!names.intree%in%mich.by.keel)
michXkeeler.tree<-drop.tip(mich.tree,to.prune)

mytree.names<-michXkeeler.tree$tip.label
michXkeeler.data<-mich.data[match(mytree.names, mich.data$name),]
namelist<-michXkeeler.data$name
namelist==mytree.names
michXkeeler.data$name== mytree.names

###########compare 2 variable models####
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<-  silv.data %>% remove_rownames %>% column_to_rownames(var="name")
keeler.data<-  keeler.data %>% remove_rownames %>% column_to_rownames(var="name")
michXkeeler.data<-  michXkeeler.data %>% remove_rownames %>% column_to_rownames(var="name")


mich1<-phyloglm(pro~pol+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=1000,full.matrix = TRUE)
summary(mich1)

silv1<-phyloglm(pro~pol+flower_time,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=1000,full.matrix = TRUE)
summary(silv1)



keeler1<-phyloglm(pro~pol+flo_time,keeler.data, keeler.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                  start.beta=NULL, start.alpha=NULL,
                  boot=10,full.matrix = TRUE)
summary(keeler1)

###WITH ONLY 2 PREDICTORS, flowering time is only significant one############but....

###ADDING TO MICHIGAN TREES MODELS####
mich2<-phyloglm(pro~pol+heigh_height+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich2)### height and flo_time are significant

mich3<-phyloglm(pro~pol+heigh_height+shade_bin+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich3)

mich4<-phyloglm(pro~pol+heigh_height+flo_time+av_fruit_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich4)
####pollination only become significant when av_fruit_time is in the model

mich5<-phyloglm(pro~pol+heigh_height+flo_time+av_fruit_time+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich5)

##########all predictors, no interactions
mich6<-phyloglm(pro~pol+heigh_height+flo_time+av_fruit_time+shade_bin+flo_type,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=10,full.matrix = TRUE)
summary(mich6)

#leave out fruit_time from model
mich7<-phyloglm(pro~pol+heigh_height+flo_time+shade_bin+flo_type,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich7) ####as i feared, pollination drops out

#So, pollination is obly significant when aveerage fruit time is in the model does this indicate an interaction?

##Does this happen in silvics?

silv2<-phyloglm(pro~pol+flower_time+av_fruit_time,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(silv2) 


s##compare to:
mich8<-phyloglm(pro~pol+av_fruit_time+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich8) ###also not



AIC(mich1)
AIC(mich2)
AIC(mich3)
AIC(mich4)####best model fruit time, flower time, height, and pollination
AIC(mich5)### 2nd best model
AIC(mich6)
AIC(mich7)
AIC(mich8)

###so pollination is significant in the best model
