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
#https://academic-oup-com.ezp-prod1.hul.harvard.edu/sysbio/article-lookup/doi/10.1093/sysbio/syp074 Garland and Ives 2010

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy<-filter(anthy, !is.na(av_fruit_time))
source("source/prune_tree.R")
is.ultrametric(pruned.by.anthy)
#############################################################################
####PGLS models for height and pollination syndrome##########################
library("phylolm")
#make $name row names
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")

#what is the data structure?
lapply(final.df, class) ### does this matter for bianary?


####full model with everything bianary#########
full.mod<-glm(pro~pol*class2*shade_bin*fruit_bin,family = binomial(link="logit"),data=final.df)
summary(full.mod)

####full  phylogentically corrected###########
full.modA<-phyloglm(pro~pol+class2+shade_bin+fruit_bin+flo_type,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(full.modA)

###full phylogenetically, with hysteranthy to include synanthy
full.modAA<-phyloglm(pro2~pol+class2+shade_bin+fruit_bin+flo_type,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot = 0, full.matrix = TRUE)
summary(full.modAA)


#########That was fun####################Nowdoit in BRMS############################################

library("brms")
library("MCMCglmm")
#https://cran.r-project.org/web/packages/brms/README.html for some guideance

#make all variable character for brms
final.df$pro<-as.character(final.df$pro)
final.df$pol<-as.character(final.df$pol)
final.df$class2<-as.character(final.df$class2)
final.df$shade_bin<-as.character(final.df$shade_bin)
final.df$fruit_bin<-as.character(final.df$fruit_bin)


##construct covarience matrix:
inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
final.df<-rownames_to_column(final.df, "name")

###Best model so far
model <- brm(pro~ pol+class2+fruit_bin+shade_bin +flo_type+ (1|name), data = final.df, 
 family = bernoulli(link="logit"), cov_ranef = list(pruned.by.anthy= A),iter=10000,
 prior = c(prior(normal(0, 5), "b"),
 prior(normal(0, 5), "Intercept"),
 prior(student_t(3, 0, 10), "sd")))
summary(model)

plot(marginal_effects(model, probs = c(0.05, 0.95)))

### not as good model
model2 <- brm(pro~ pol+class2 +fruit_bin+shade_bin+ (1|name), data = final.df, 
                    family = bernoulli(link="logit"), cov_ranef = list(pruned.by.anthy= A),iter=10000,
                    prior = c(prior(cauchy(0, 5), "b"),
                              prior(cauchy(0, 5), "Intercept"),
                              prior(student_t(3,0, 10), "sd")))   


summary(model2)

library(shinystan)                   
launch_shiny(model, rstudio = getOption("shinystan.rstudio"))


#####compute phylo.d########################################################
#final.df<-rownames_to_column(final.df  ,var="rowname")
#phylo.d(data = final.df,phy = pruned.by.anthy, names.col = name, binvar = pro, permut = 1000, rnd.bias = NULL)

myy<-final.df$pro


