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


#####compute phylo.d########################################################
#final.df<-select_(final.df,"name","pro","pol")
#final.df<-filter(final.df, name== c(namelist2))
#which(final.df$name%in%pruned.by.anthy$tip.label)
#which(pruned.by.anthy$tip.label%in%final.df$name)
#pruned.by.anthy$node.label<-""
#phylo.d(data = final.df,phy = pruned.by.anthy, names.col = tip.label, binvar = pro, permut = 1000, rnd.bias = NULL)
#d<-comparative.data(pruned.by.anthy,final.df,name, na.omit = FALSE)
#d

#phylo.d(data = final.df,phy = pruned.by.anthy, names.col = name, binvar = pro, permut = 1000, rnd.bias = NULL)
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

####add flower type
mod4<-glm(pro~pol+class2+flower_class,family = binomial(link="logit"),data=final.df)
summary(mod4)
inv.logit(coef(mod4)[1])#= 0.06293043  if you are dioecious insect shrub
inv.logit(coef(mod4)[1]+coef(mod4)[2])#=0.7539598 if you are dioecious wind shrub
inv.logit(coef(mod4)[1]+coef(mod4)[3])#= 0.1259112 if you are dioecious insect tree
inv.logit(coef(mod4)[1]+coef(mod4)[4])#=0.002195789 if you are a monecious insect shrub
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[3])#= 0.8679515 if you are a dioecious wind tree
inv.logit(coef(mod4)[1]+coef(mod4)[2]+coef(mod4)[3]+coef(mod4)[4])#=0.1772165 if you are a monecious wind tree


#mod4a<-phyloglm(pro~pol+class2+flower_class,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                #start.beta=NULL, start.alpha=NULL,
                #boot = 0, full.matrix = TRUE)
#summary(mod4a)



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
help(brm)
model_simple <- brm(pro ~ pol +class2+ (1|name), data = final.df, 
                    family = bernoulli(link = "logit"), cov_ranef = list(pruned.by.anthy = A), iter= 10000)
###priors?
#cov_ranef = list(pruned.by.anthy = A),
pairs(model_simple)
summary(model_simple)
library(shinystan)                   
launch_shiny(model_simple, rstudio = getOption("shinystan.rstudio"))

