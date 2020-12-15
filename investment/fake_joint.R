## Started 15 Dec 2020
## Dan's fake data and model buidling for a joing model about FLSs in Prunus


# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(rstan)

set.seed(40)
setwd("~/Documents/git/ospree/analyses/ranges")

#--------------------------------------
# Set up the FLS model - 
# FLS ~  agrand a[sp] + doy + sigma_y (doy is a covvariate controling for temporality of samplinh)

# Parameters
a_grand <- 9
sigma_y <- 1

beta_doy<-0.025

n <- 200 # number of replicates per sp x study (may eventually want to draw this from a distribution to make data more realistic)
nsp <- 13 # number of species

### This will be a fixed effects model but I think we need some mua_sp to create some variation around our species estimates
sigma_asp <- 0.7
mua_sp <- rnorm(nsp, 0, sigma_asp)

# Set up the data ...
simFLS <- data.frame(sp=rep(1:nsp, each=200), mua_sp=rep(mua_sp, each=200))
simFLS$doy<-rnorm(nrow(simFLS), 120,20)

simFLS$FLS <- a_grand + simFLS$mua_sp+beta_doy*simFLS$doy + rnorm(nrow(simFLS), 0, sigma_y)
hist(simFLS$FLS)


ggpubr::ggboxplot(simFLS,'sp','FLS')

library(lme4)
# Fixed effects should be very close to coefficients used in simulating data
summary(lme1 <- lmer(FLS ~ doy+(1|sp), data = simFLS)) 
plot(ranef(lme1)$sp[,], mua_sp)

#--------------------------------------
# Now simulate the pdsi side
# pdsi ~ FLS_[sp] not quite sure ghow to do this



