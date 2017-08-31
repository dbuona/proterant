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

####read in tree
my.tree<-read.tree( "pruned_for_mich.tre")
my.data<-read.csv( "mich_data_full.csv", header=TRUE)

my.data<-dplyr::select(my.data, name, pro,pol,heigh_height, shade_bin, av_fruit_time, flo_time)
colnames(my.data)<- c("name","y","c","d","e","f","g")
write.csv(my.data, "data_example.csv")

###diguise variable


full.mod<-glm(y~c+d+e+f+g,family = binomial(link="logit"),data=my.data)
summary(full.mod)

my.data<-  my.data %>% remove_rownames %>% column_to_rownames(var="name")
full.modB<-phyloglm(y~c+d+e+f+g,my.data, my.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 100,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)
summary(full.modB)
####pollination drops out without a fruit_time metric. Why

###Bayesian
library("brms")
library("MCMCglmm")

my.data<-rownames_to_column(my.data, "name")
inv.phylo <- MCMCglmm::inverseA(my.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###bayesian and continuous-- main model###############
modelcont <- brm(y~c+d+e+f+g+ (1|name), data = my.data, 
                 family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
                 prior = c(prior(normal(0, 5), "b"),
                           prior(normal(0, 5), "Intercept"),
                           prior(student_t(3, 0, 5), "sd"))) 

summary(modelcont)



