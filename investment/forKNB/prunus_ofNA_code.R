##All American Prunus  4KNB

###housekeeping
rm(list=ls()) ## remove everything in R memory
options(stringsAsFactors = FALSE) # strings are not factors
options(mc.cores = parallel::detectCores()) ## parallelize chains for stan/brms

###install/load packages
library(dplyr)
library(ggplot2)
library(rstan)
library(brms)
library(phytools)
library(ape)
library(tidybayes)
library(ggpubr)

setwd("") ## set your working directory

d<-read.csv("prunus_final.csv") ##read in data
tr<-read.tree("prunus_tree.tre") ## read in tree

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)} ## fucntion for zscoring

d$meanpdsi.z<-zscore(d$meanpdsi) ## zscore predictors
d$inflor.z<-zscore(d$inflor_high)


A <- ape::vcv.phylo(tr) ### set up VCV matrix with phylogenetic distances

###run model in brms
FNAordz.phylo<-brm(FLSnum~inflor.z*meanpdsi.z +(1|specificEpithet)+(1|gr(species, cov = A)),
                    data=d,
                    family=cumulative("logit"),
                    data2 = list(A = A),control=list(adapt_delta=0.99),
                    warmup=6000,iter=8000)


###gather the effect size estimates
output<-FNAordz.phylo %>%
  spread_draws(b_inflor.z ,b_meanpdsi.z,`b_inflor.z:meanpdsi.z`)
colnames(output)
output <-output %>% tidyr::gather("var","estimate",4:6)
#####

###Plot 4a#############################
ggplot(output,aes(y = var, x = estimate)) +
  stat_pointinterval(.width=c(.5,.89))+ggthemes::theme_few()+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_y_discrete(limits = c("b_inflor.z:meanpdsi.z","b_meanpdsi.z","b_inflor.z"),labels=c("PDSI X infloresence size","PDSI","inflorescence size"))+
  ylab("")+xlab("standardized effect size estimate")

## plot 4b with marginal effects
p1<-plot(conditional_effects(FNAordz.phylo, "meanpdsi.z", ordinal = TRUE,prob = .5,plot=FALSE))
p3<-plot(conditional_effects(FNAordz.phylo, "inflor.z", ordinal = TRUE,prob=.5,plot=FALSE))

##make prettier and label
p1<-p1[[1]]+ggthemes::theme_few()+scale_y_discrete(name="FLS",labels=c("flowers after leaves","flowers with leaves","flowers before/with leaves","flowers before leaves"))+xlab("PDSI")
p3<-p3[[1]]+ggthemes::theme_few()+ylab("")+xlab("inflorescence size")+theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())

### plot side by side
ggpubr::ggarrange(p1,p3,common.legend = TRUE,ncol=2,legend="bottom",widths = c(.8,.5))
