
# Working on Bayes PGLS, following along here #
# https://github.com/Auerilas/BayesPGLS/blob/master/PGLS.py #
# Same as Will Pearse's code I think ... ##

## Updated on 17 May 2016 (Happy Birthday! and Happy J Tenure day!) ##
## By Lizzie ##

## This file now adapts the python file and seems to run ##
## I made a copy of this file (with minor directory adjustments) in .. ##
# teaching/gelmanhill/BayesPGLS ##

library(phytools)
library(rstan)
library("shinystan")

setwd("~/Documents/git/teaching/stan/phylogeny")
shorebirdTraits <- read.csv('input/shorebirdData.csv', header=TRUE, row.names=1)
shorebirdVCVdat <- read.csv('input/shorebirdVCV.csv', skip=1, header=FALSE)
shorebirdVCV <- as.matrix(shorebirdVCVdat)

Lmat <- matrix(rep(1), nrow = nrow(shorebirdVCV), ncol = ncol(shorebirdVCV))
diag(Lmat) <- 0

stdX <- scale(shorebirdTraits$M.Mass, center=TRUE, scale=TRUE)
stdY <- scale(shorebirdTraits$Egg.Mass, center=TRUE, scale=TRUE)

## build up data for stan model
N <- nrow(shorebirdTraits)
X <- as.vector(stdX)
K <- 1
V <- shorebirdVCV
y <- as.vector(stdY)

fit.pgls <- stan("stan/pgls_lemoine.stan", data=c("N","X", "K", "Lmat", "V", "y"), iter=2000, chains=4)
launch_shinystan(fit.pgls)

## compare with caper
library(caper)
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
caper.mod <- pgls(scale(Egg.Mass) ~ scale(M.Mass), data=shorebird, lambda='ML')

summary(fit.pgls)
