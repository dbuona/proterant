rm(list=ls()) 
options(stringsAsFactors = FALSE)
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
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
source("source/prune_tree.R")

#####compute phylo.d########################################################
#final.df<-select_(final.df,"name","pro","pol")
#final.df<-filter(final.df, name== c(namelist2))
#which(final.df$name%in%pruned.by.anthy$tip.label)
#which(pruned.by.anthy$tip.label%in%final.df$name)
#pruned.by.anthy$node.label<-""
#phylo.d(data = final.df,phy = pruned.by.anthy, names.col = tip.label, binvar = pro, permut = 1000, rnd.bias = NULL)
#d<-comparative.data(pruned.by.anthy,final.df,rownames, na.omit = FALSE)
#d
#phylo.d(d,binvar=pro)

#############################################################################
####PGLS models for height and pollination syndrome##########################
library("phylolm")
#make $name row names
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")

mod1<-glm(pro~pol+heigh_height,family = binomial(link="logit"),data=final.df)
summary(mod1)
inv.logit(-2.00021) #=0.1191809
inv.logit(-2.0002+1.35708) #=0.3445461
###but this is only if height is zero....still need help to interpret

inv.logit(coef(mod1)[1]+coef(mod1)[3]*mean(final.df$heigh_height))##gellman does this

mod1a<-phyloglm(pro~pol+height10,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 30, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod1a)
######PGLS models for height as discrete pollin syndrome
mod2<-glm(pro~pol+class,family = binomial(link="logit"),data=final.df)
summary(mod2)
inv.logit(coef(mod2)[1])#= 1902566 % chance if you are a insect shrub
inv.logit(coef(mod2)[1]+coef(mod2)[2])#=0.5395037 if you are wind shrub
inv.logit(coef(mod2)[1]+coef(mod2)[3])#=0.1725329 if you are a insect small tree
inv.logit(coef(mod2)[1]+coef(mod2)[2]+coef(mod2)[3])#=0.5095634 if you are wind small tree

inv.logit(coef(mod2)[1]+coef(mod2)[2]+coef(mod2)[3]+coef(mod2)[3])#= 0.4795543if you are a large wind pollen

mod2a<-phyloglm(pro~pol+class,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod2a)
inv.logit(coef(mod2a)[1])#= 0.2493044 if you are insect shrub
inv.logit(coef(mod2a)[1]+coef(mod2a)[2])#=0.4614641 if you are wind shrub
inv.logit(coef(mod2a)[1]+coef(mod2a)[3])#=0.2040091  if you are a insect small tree
inv.logit(coef(mod2a)[1]+coef(mod2a)[2]+coef(mod2a)[3])#=0.3980618 if you are wind small tree
##how to do large tree?
inv.logit(coef(mod2a)[1]+coef(mod2a)[2]+coef(mod2a)[3]+coef(mod2a)[3])#=0.3379048 if you are a large wind pollen tree


##same as above for phylogenetically adjusted

######PGLS models for shrub/tree bianary and pollin syndrome
mod3<-glm(pro~pol+class2,family = binomial(link="logit"),data=final.df)
summary(mod3)
inv.logit(coef(mod3)[1])#= 0.1518854 if you are insect shrub
inv.logit(coef(mod3)[1]+coef(mod3)[2])#=0.4228354 if you are wind shrub
inv.logit(coef(mod3)[1]+coef(mod3)[3])#= 0.1946361 if you are insect tree
inv.logit(coef(mod3)[1]+coef(mod3)[2]+coef(mod3)[3])#= 0.4971456 if you are wind tree

mod3a<-phyloglm(pro~pol+class2,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod3a)
inv.logit(coef(mod3a)[1])#= 0.1657828  if you are insect shrub
inv.logit(coef(mod3a)[1]+coef(mod3a)[2])#=0.2775916 if you are wind shrub
inv.logit(coef(mod3a)[1]+coef(mod3a)[3])#= 0.2203852  if you are insect tree
inv.logit(coef(mod3a)[1]+coef(mod3a)[2]+coef(mod3a)[3])#= 0.3534184  if you are wind tree


#########That was fun####################Nowdoit in BRMS####
#install.packages("brms") #Error in install.packages : cannot remove prior installation of package ‘scales’ ## try this again on different computer
library("brms")
library("MCMCglmm")


model_simple <- brm(pro ~ pol +class2+ (1|pruned.by.anthy), data = final.df, 
                    family = bernoulli(link = "logit"), cov_ranef = list(phylo = A))

                   


