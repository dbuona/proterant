### put all you models in one place
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library("rstan")

library(phytools)
library(ape)
library(lubridate)
library(stringr)
library("tidybayes")
library(raster)
library(phylolm)
library(tibble)

require(mapdata); require(maptools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/Documents/git/proterant/investment/Input")

##read in cleaned data
d.flo<-read.csv("input_clean/FLS_clean.csv")
d.pdsi<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi.csv")
d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")
d.phen<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruit_phen.csv")
d.cold<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi_wint.csv")

dater.flo<-d.flo %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanFLS=mean(bbch.short),sdFLS=sd(bbch.short))
dater.pheno<-d.flo %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanFLO=mean(doy),sdFLO=sd(doy))
dater.pdsi<-d.pdsi %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanPDSI=mean(pdsi,na.rm=TRUE),sdPDSI=sd(pdsi,na.rm=TRUE),minPDSI=mean(pdsi.min,na.rm=TRUE),sdmin=sd(pdsi.min,na.rm=TRUE))
dater.petal<-d.petal %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanPETAL=mean(pental_lengh_mm,na.rm=TRUE),sdPETAL=sd(pental_lengh_mm,na.rm=TRUE))
d.fruit<-filter(d.fruit,fruit_type=="fleshy")
dater.fruit<-d.fruit %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanFRUIT=mean(fruit_diam_mm,na.rm=TRUE),sdFRUIT=sd(fruit_diam_mm,na.rm=TRUE))
dater<-left_join(dater.flo,dater.pdsi)
dater<-left_join(dater,dater.petal)
dater<-left_join(dater,dater.fruit)
dater<-left_join(dater,dater.pheno)

tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre")
is.ultrametric(tree)
tree<-compute.brlen(tree, method = "Grafen")
## make ultrametric

names.intree<-tree$tip.label # names the names
namelist<-unique(d.flo$specificEpithet)

to.prune<-which(!names.intree%in%namelist) #prun the tree
pruned.by<-drop.tip(tree,to.prune)
plotTree(pruned.by)# this is the tree
mytree.names<-pruned.by$tip.label

final.df<-dater[match(mytree.names, dater$specificEpithet),]
namelist2<-final.df$specificEpithet
namelist2==mytree.names
final.df$specificEpithet== mytree.names
final.df<-  final.df %>% remove_rownames() %>% column_to_rownames(var="specificEpithet")

pruno.phylo.lm<-phylolm(meanPETAL~meanFLS*meanFLO+meanPDSI,data=final.df,phy=pruned.by)
pruno.lm<-lm(meanFLS~meanPDSI,data=final.df)


summary(pruno.phylo.lm)
summary(pruno.lm)


###mean model############
mod1<-brm(meanFLS~meanPDSI,data=dater)
A <- ape::vcv.phylo(pruned.by)

mod1phylo<-brm(meanFLS~+meanPDSI + (1|gr(specificEpithet, cov = A)),
                   data=dater,
                   data2 = list(A = A),
                   warmup=3000,iter=4000,control=list(adapt_delta=.995))


fixef(mod1,probs = c(.25,.75))
fixef(mod1phylo,probs = c(.25,.75))
fixef(mod1phylo.mi,probs = c(.25,.75))

mod1a<-brm(meanFLS|mi(sdFLS)~meanPDSI,data=dater)

mod1phylo.mi<-brm(meanFLS|mi(sdFLS)~meanPDSI + (1|gr(specificEpithet, cov = A)),
               data=dater,
               data2 = list(A = A),
               warmup=3000,iter=4000,control=list(adapt_delta=.995))


mod2<-brm(meanFLS| mi(sdFLS) ~ me(meanPDSI,sdPDSI),data=dater, control=list(adapt_delta=.95))

#####measurement error)
