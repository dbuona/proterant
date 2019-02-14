# Working on Bayes PGLS, following along here #
# https://github.com/Auerilas/BayesPGLS/blob/master/PGLS.py #
# Same as Will Pearse's code I think ... ##

## Updated on 17 May 2016 (Happy Birthday! and Happy J Tenure day!) ##
## By Lizzie ##

## This file now adapts the python file and seems to run ##
## I made a copy of this file (with minor directory adjustments) in .. ##
# teaching/gelmanhill/BayesPGLS ##


#### Dan B began modifying this on Feb 11 2019

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")

library(phytools)
library(rstan)
library("shinystan")
library(ape)
library(tidyverse)


mich.data<-read.csv("datasheets_derived/MTSV_USFS/mich_data_full_clean.csv")
mich.tree<-read.tree("datasheets_derived/MTSV_USFS/pruned_for_mich.tre")
trait.data<-read.csv("USDA_traitfor_MTSV.csv",header=TRUE)

###clean data
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)
###make the tree work 
mich.tree$node.label<-NULL
mich.data<-left_join(mich.data,trait.data,by="name")

######prune the tree for drought modeling
mich.data<-filter(mich.data,!is.na(min._precip))

####prune tree to match reduced dataset
names.intree<-mich.tree$tip.label
namelist<-unique(mich.data$name)
to.prune<-which(!names.intree%in%namelist)
mich.tree.droughtprune<-drop.tip(mich.tree,to.prune)
mytree.names<-mich.tree.droughtprune$tip.label
### the smaller new tree is called mich.tree.droughtprune

####recenter here after you prune the list
###Rescale predictors this makes it so you can compare binary to continous data
#mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
#mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))
#mich.data$dev_time_cent<-(mich.data$dev.time-mean(mich.data$dev.time))/(2*sd(mich.data$dev.time))
#mich.data$tol_cent<-(mich.data$shade_bin-mean(mich.data$shade_bin))/(2*sd(mich.data$shade_bin))
mich.data$precip_cent<-(mich.data$min._precip-mean(mich.data$min._precip))/(2*sd(mich.data$min._precip))

####now you have you data a you matrix

michVCV <- vcv.phylo(mich.tree.droughtprune, cor=TRUE) #this could extract my vVC. from the file

Lmat <- matrix(rep(1), nrow = nrow(michVCV), ncol = ncol(michVCV))
diag(Lmat) <- 0

#dont neet this
#stdX <- scale(shorebirdTraits$M.Mass, center=TRUE, scale=TRUE)
#stdY <- scale(shorebirdTraits$Egg.Mass, center=TRUE, scale=TRUE)

## build up data for stan model
N <- nrow(mich.data)
flo_time <- as.vector(mich.data$flo_time)
K <- 1
V <- michVCV
y <- as.vector(mich.data$pro2)

fit.pgls <- stan("..//stan/pgls_hysteranthy.stan", data=c("N","flo_time", "K", "Lmat", "V", "y"), iter=2000, chains=4)
launch_shinystan(fit.pgls)

## compare with caper
library(caper)
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
caper.mod <- pgls(scale(Egg.Mass) ~ scale(M.Mass), data=shorebird, lambda='ML')

summary(fit.pgls)
