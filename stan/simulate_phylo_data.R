
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()
library(phylolm)
library(extraDistr)
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
library(car)
library(ggstance)
library(broom)
library(brms)
library(rstan)

#library(devtools)
#install_github("rmcelreath/rethinking")
library(rstan)
setwd("~/Documents/git/proterant/Input/datasheets_derived/MTSV_USFS/")


set.seed(613)



tre <-read.tree("pruned_for_mich.tre")
tre$node.label<-NULL
lab<-tre$tip.label 
### simulates 2 continues variables with differnt phylogenetic clustering themselves

?rTrait()
x = rTrait(n=2,phy=tre,model="lambda")
y = rbinTrait(n=1, phy=tre, beta=c(-.5,1.8), alpha=0, X=x)
dat<-as.data.frame(cbind(x,y))


set.seed(123456)
tre = rtree(50)
x = rTrait(n=1,phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=0 ,X=X)
dat = data.frame(trait01 = y, predictor = x)
fit = phyloglm(trait01~predictor,phy=tre,data=dat,boot=100)
summary(fit)

?phyloglm()
summary(mod1)
dat2<-rownames_to_column(dat,var = "name")

inv.phylo <- MCMCglmm::inverseA(tre, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

###continuous models
brms.try <- brm(y~V1+V2+(1|name), data = dat2, 
                             family = gaussian(), cov_ranef = list(name= A),iter=3000) 

summary(brms.try)
