### Started 1 April 2019 - Cat
## Building stan models to assess impact of false springs

### Main model:
#meri.mod <- brm(meristem ~ tx*chill1 + tx*chill2 + (tx*chill1 + tx*chill2 | species),
 #               data=chill.stan, family=binomial(link="logit"), iter=4000, warmup=2500, 
  #              control=list(max_treedepth=15, adapt_delta=0.99))
#save(meri.mod, file="~/Documents/git/chillfreeze/analyses/stan/meristem_brms.Rdata")


# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Load libraries
library(RColorBrewer)
library(rstan)
library(dplyr)
library(broom)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory
setwd("~/Documents/git/chillfreeze/analyses")

## load the model
load("stan/meristem_brms.Rdata")

chill.stan <- read.csv("output/clean_dvr_traits.csv") # load the data
chill.stan <- chill.stan[!is.na(chill.stan$meristem),]

chill.stan$species.name <- NA
chill.stan$species.name <- ifelse(chill.stan$species=="ACESAC", "Acer saccharinum", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="ALNRUG", "Alnus rugosa", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="BETPAP", "Betula papyrifera", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="BETPOP", "Betula populifolia", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="CORRAC", "Cornus racemosa", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="SALPUR", "Salix purpurea", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="SORAME", "Sorbus americana", chill.stan$species.name)
chill.stan$species.name <- ifelse(chill.stan$species=="VIBDEN", "Viburnum dentatum", chill.stan$species.name)


#### Now for mu plots based of bb_analysis/models_stan_plotting.R ###
figpath <- "figures"
figpathmore <- "meristem_brms" ### change based on model

source("exp_muplot_brms.R")
cols <- adjustcolor("indianred3", alpha.f = 0.3) 
my.pal <- rep(brewer.pal(n = 10, name = "Paired"), 8)
# display.brewer.all()
alphahere = 0.4

xlab <- "Model estimate of change in shoot apical meristem damage"

spp <- unique(chill.stan$species)

modelhere <- meri.mod

tx <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 2] %>% ## here we make find the posterior distributions and means for each predictor
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select(mean, `25%`, `75%`) ### can change according to uncertainty intervals you want
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("tx", "[", i, "]", sep="")
}
tx$parameter<-new.names
chill1 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 3] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("chill1", "[", i, "]", sep="")
}
chill1$parameter<-new.names
mod.ranef<-full_join(tx, chill1)
chill2 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 4] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("chill2", "[", i, "]", sep="")
}
chill2$parameter<-new.names
mod.ranef <- full_join(mod.ranef, chill2)
txchill1 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 5] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("tx:chill1", "[", i, "]", sep="")
}
txchill1$parameter<-new.names
mod.ranef<-full_join(mod.ranef, txchill1)
txchill2 <- coef(modelhere, prob=c(0.25, 0.75))$species[, c(1, 3:4), 6] %>%
  as.data.frame() %>%
  round(digits = 2) %>% 
  rename(mean = Estimate) %>%
  rename(`25%` = Q25) %>%
  rename(`75%` = Q75) %>%
  dplyr::select( mean, `25%`, `75%`) 
new.names<-NULL
for(i in 1:length(spp)){
  new.names[i]<-paste("tx:chill2", "[", i, "]", sep="")
}
txchill2$parameter<-new.names
mod.ranef<-full_join(mod.ranef, txchill2)

modoutput <- tidy(modelhere, prob=c(0.5))

muplotfx(modelhere, "", 8, 8, c(0,5), c(-8, 12) , 12.5, 3.5)

