####Fake data for hysteranty
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
library(msm) 
set.seed(73)
##simulate tree

tre = rtree(200)
tre$tip.label
tre<-chronoMPL(tre)
is.ultrametric(tre)
###data
ntot = 200

pro<-rbinom(ntot, 1, 0.5)
pol<-rbinom(ntot, 1, 0.5)
tol<-rbinom(ntot, 1, 0.5)
height<-rtnorm(ntot, 20, 5, lower=3, upper=40)
flower<-rtnorm(ntot, 5, 2, lower=2, upper=12)

mm <-as.data.frame(model.matrix(~(pro+pol+tol+height+flower), data.frame(pro+pol+tol+height+flower)))
mm$name<-tre$tip.label
mm<-dplyr::select(mm,name,pro,pol,tol,height,flower)

data<-mm[match(tre$tip.label, mm$name),]
mytree.names<-tre$tip.label
namelist2<-data$name
namelist2==mytree.names
data$name== mytree.names
data<-  data %>% remove_rownames %>% column_to_rownames(var="name")

mod<-glm(pro~pol+tol+height+flower,family=binomial(link="logit"), data=data)
summary(mod)   
    
full.modB<-phyloglm(pro~pol+tol+height+flower, data,tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)
summary(full.modB)
####

library("brms")
library("MCMCglmm")


data<-rownames_to_column(data, "name")
inv.phylo <- MCMCglmm::inverseA(tre, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###bayesian and continuous-- main model###############
modelcont <- brm(pro~ pol+tol+height+flower+(1|name), data = data, 
                 family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
                 prior = c(prior(normal(0, 5), "b"),
                           prior(normal(0, 5), "Intercept"),
                           prior(student_t(3, 0, 5), "sd"))) 


 