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


###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
source("source/prune_tree.R")

#############################################################################
####PGLS models for height and pollination syndrome##########################
library("phylolm")
#make $name row names
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")

######PGLS models for shrub/tree bianary and pollin syndrome (this is model i think i care about the most #################
mod3<-glm(pro~pol+class2,family = binomial(link="logit"),data=final.df)
summary(mod3)
inv.logit(coef(mod3)[1])
inv.logit(coef(mod3)[1]+coef(mod3)[2])
inv.logit(coef(mod3)[1]+coef(mod3)[3])
inv.logit(coef(mod3)[1]+coef(mod3)[2]+coef(mod3)[3])

mod3a<-phyloglm(pro~pol+class2,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod3a)
inv.logit(coef(mod3a)[1])#if you are insect shrub
inv.logit(coef(mod3a)[1]+coef(mod3a)[2])
inv.logit(coef(mod3a)[1]+coef(mod3a)[3])
inv.logit(coef(mod3a)[1]+coef(mod3a)[2]+coef(mod3a)[3])


#########################################################add flower type
mod4<-glm(pro~pol+class2+flo_type,family = binomial(link="logit"),data=final.df)
summary(mod4)
inv.logit(coef(mod4)[1])#  if you are perfect insect shrub
inv.logit(coef(mod4)[1]+coef(mod4)[2]) #perfect wind shrub
inv.logit(coef(mod4)[1]+coef(mod4)[3]) #perfect insect tree
inv.logit(coef(mod4)[1]+coef(mod4)[4])#  monecious insect  shrub
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[3]) #perfect wind tree
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[3]+coef(mod4)[4]) # mono, wind, tree
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[4]) # mono wind shrub
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[3]+coef(mod4)[5]) #dio wind tree
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[5]) #dio wind shrub


mod4a<-phyloglm(pro~pol+class2+flo_type,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod4a)
inv.logit(coef(mod4a)[1])#  if you are perfect insect shrub
inv.logit(coef(mod4a)[1]+coef(mod4a)[2]) #perfect wind shrub
inv.logit(coef(mod4a)[1]+coef(mod4a)[3]) #perfect insect tree
inv.logit(coef(mod4a)[1]+coef(mod4a)[4])#  monecious insect  shrub
inv.logit(coef(mod4a)[1]+coef(mod4a)[2]+coef(mod4a)[3]) #perfect wind tree
inv.logit(coef(mod4a)[1]+coef(mod4a)[2]+coef(mod4a)[3]+coef(mod4a)[4]) # mono, wind, tree
inv.logit(coef(mod4a)[1]+coef(mod4a)[2]+coef(mod4a)[4]) # mono wind shrub
inv.logit(coef(mod4a)[1]+coef(mod4a)[2]+coef(mod4a)[3]+coef(mod4a)[5]) #dio wind tree
inv.logit(coef(mod4a)[1]+coef(mod4)[2]+coef(mod4a)[5])



##metrics
#wind pollinated shrubs
dim(filter(final.df, pol==1 & class=="0shrub"))
dim(filter(final.df, pol==1 & class2==1))
dim(filter(final.df,pol==1 &flo_type=="0perfect"& flo_type=="0perfect")) #just elms
dim(filter(final.df,pol==1 &flo_type=="0perfect" & class=="0shrub")) ##doesnt exist 


#stop("just stop for now")
#########That was fun####################Nowdoit in BRMS####
library("brms")
library("MCMCglmm")

pruned.by.anthy$node.label<-""
###make ultrametric (using mean path length smoothing, could also try penalized maximum likelihood with chronos())
is.ultrametric(pruned.by.anthy)
plot(pruned.by.anthy)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)
plot(pruned.by.anthy)
##construct covarience matrix
inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
final.df<-rownames_to_column(final.df, "name")

model_simple <- brm(pro~ pol+class2 + (1|name), data = final.df, 
 family = bernoulli(link="logit"), cov_ranef = list(pruned.by.anthy= A),
 prior = c(prior(normal(0, 5), "b"),
 prior(normal(0, 5), "Intercept"),
 prior(student_t(3, 0, 10), "sd")))                    
######need help determining priors--and interpretting output. This provides very different estimate than plgs
pairs(model_simple)
summary(model_simple)

library(shinystan)                   
launch_shiny(model_simple, rstudio = getOption("shinystan.rstudio"))

#####compute phylo.d########################################################
#final.df<-rownames_to_column(final.df  ,var="rowname")
#phylo.d(data = final.df,phy = pruned.by.anthy, names.col = rowname, binvar = pro, permut = 1000, rnd.bias = NULL)


