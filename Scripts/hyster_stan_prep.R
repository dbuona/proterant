###Dan's attempt at getting his hysternathy data into stan
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(ape)
library(phytools)
library(gbm)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(tidyverse)


#https://academic-oup-com.ezp-prod1.hul.harvard.edu/sysbio/article-lookup/doi/10.1093/sysbio/syp074 Garland and Ives 2010

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy<-filter(anthy, !is.na(av_fruit_time))
source("source/prune_tree.R") ###this prunes the tree and cleans the data
###Tree is pruned.by.anthy
###data is final.df

##### make VCV matrix based on tree
#library("MCMCglmm")

covv<-vcv.phylo(pruned.by.anthy, model = "Brownian", cor = FALSE)

#inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
#A <- solve(inv.phylo$Ainv)
#rownames(A) <- rownames(inv.phylo$Ainv)

## prepare data (square matrix of 1s and 0s in main diagonal)
Lmat <- matrix(rep(1), nrow = nrow(covv), ncol = ncol(covv))
diag(Lmat) <- 0

## build up data for stan model
N <- nrow(final.df) # nspecies (analysis units)
pol<- as.vector(final.df$pol) # predictor
class2<- as.vector(final.df$class2) # predictor
shade_bin <- as.vector(final.df$shade_bin) # predictor
fruit_bin<-as.vector(final.df$fruit_bin) # predictor
flo_type<-as.vector(final.df$flo_type) #predictor
sp=final.df$name ###not sure if this is right
n_sp=length(unique(final.df$name))

K <- 5 #number of predictors?
V <- as.matrix(covv) # variance covariance matrix
y <- as.vector(final.df$pro) # response variable
library(rstan)
fit.pgls <- stan("hyster.stan", data=c("N","pol","class2","shade_bin","fruit_bin","flo_type", "K", "Lmat", "V", "y"), iter=2000, chains=4)
launch_shinystan(fit.pgls)
