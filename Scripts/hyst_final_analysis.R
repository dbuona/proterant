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

#####add a new column for a adjusting for red acorn time
mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
mich.data$fruiting[mich.data$fruiting==19]<-7

silv.data$fruiting<-NA
silv.data$fruiting<-silv.data$av_fruit_time
silv.data$fruiting[silv.data$fruiting==21]<-9

###########compare 2 variable models####
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<-  silv.data %>% remove_rownames %>% column_to_rownames(var="name")
keeler.data<-  keeler.data %>% remove_rownames %>% column_to_rownames(var="name")
michXkeeler.data<-  michXkeeler.data %>% remove_rownames %>% column_to_rownames(var="name")

mich.2<-phyloglm(pro~pol+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich2)

###ADDING TO MICHIGAN TREES MODELS####
mich3<-phyloglm(pro~pol+heigh_height+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich3)### height and flo_time are significant

mich3a<-phyloglm(pro~pol+flo_time+fruiting,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich3a)

cor(mich.data$fruiting,mich.data$flo_time)

mich4<-phyloglm(pro~pol+heigh_height+shade_bin+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich4)

mich4a<-phyloglm(pro~pol+heigh_height+flo_time+fruiting,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(mich4a)
####pollination only become significant when av_fruit_time is in the model

mich5<-phyloglm(pro~pol+heigh_height+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5)


################Other datasets:

silv2<-phyloglm(pro~pol+flower_time,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(silv2)


silv3<-phyloglm(pro~pol+flower_time+av_fruit_time,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(silv3) 

silv3a<-phyloglm(pro~pol+flower_time+fruiting+heigh_height,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)
summary(silv3a)

##Cant estimate it with height

keeler2<-phyloglm(pro~pol+flo_time,keeler.data, keeler.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                  start.beta=NULL, start.alpha=NULL,
                  boot=10,full.matrix = TRUE)
summary(keeler2)

michXkeeler2<-phyloglm(pro~pol+flo_time,michXkeeler.data, michXkeeler.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                  start.beta=NULL, start.alpha=NULL,
                  boot=10,full.matrix = TRUE)
summary(michXkeeler2)
#### Can't really sink the models, I'll maybe try it with texas tomorrow, in the meantime

#plotting.
###the tree and variable:
#plotting

summary(mich5)

coef(mich5)
est<-as.data.frame(coef(mich5))
est<-rownames_to_column(est, "name")
ints<-as.data.frame(confint(mich5,level = 0.95))
ints<-rownames_to_column(ints, "name")
colnames(ints)[2] <- "low"
colnames(ints)[3] <- "high"
colnames(est)[2] <- "estimate"
foo<-left_join(est,ints)
foo<-filter(foo,estimate<10)
ggplot(foo,aes(estimate,name))+geom_point()+geom_segment(aes(y=name,yend=name,x=low,xend=high))+ggtitle("Main effects of predictors on Hysteranthy")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")

boot<-read.csv("mich5bootoutput.csv",header=TRUE)
head(boot)
boot<-dplyr::filter(boot,Coefficients!="(Intercept)")
ggplot(boot,aes(Estimate,Coefficients))+geom_point()+geom_segment(aes(y=Coefficients,yend=Coefficients,x=lowerbootCI,xend=upperbootCI))+ggtitle("Main effects of predictors on Hysteranthy with bootstrap")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")
