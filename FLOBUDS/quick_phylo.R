####main cleaning file for MTSV hysteranthy dataset

rm(list=ls()) 
options(stringsAsFactors = FALSE)
#graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS")
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(boot)
library(phylolm)
library(dplyr)


treee<-read.tree("..//Data/Vascular_Plants_rooted.dated.tre") #read in tree
names.intree<-treee$tip.label # names the names


namelist<-c("Prunus_pensylvanica","Prunus_virginiana","Acer_rubrum","Acer_pensylvanicum",
          "Ilex_mucronata","Ilex_verticillata","Vaccinium_corymbosum","Corylus_cornuta",
          "Comptonia_peregrina") # these are species I want

to.prune<-which(!names.intree%in%namelist) #prun the tree
pruned.by.anthy<-drop.tip(treee,to.prune)
plot(pruned.by.anthy)# this is the tree

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label # did i get them all

intersect(namelist,mytree.names) #yes


dat<-read.csv("flobudsdata.use.csv",header = TRUE) ## reaadin in pheno data
dat<-dplyr::filter(dat,Chill==1) # subset to like conditions High chill
dat<-dplyr::filter(dat,Force=="W") # high force



dat$FLS<- dat$flo_day-dat$budburst.9. # calculate FLSs

d<- dat %>% dplyr::group_by(GEN.SPA) %>% dplyr::summarise(meanFLS=mean(FLS,na.rm=TRUE)) ## take species mean

d<-dplyr::filter(d,!GEN.SPA%in% c("ACE.SAC","BET.ALL","VIB.ACE","AME.SPP")) # remvoe species with no FLS data

d$taxa<-c("Acer_pensylvanicum","Acer_rubrum","Comptonia_peregrina","Corylus_cornuta","Ilex_mucronata","Ilex_verticillata",
          "Prunus_pensylvanica","Prunus_virginiana","Vaccinium_corymbosum") # add a cloumn matching the tree name formate

### put it in same order as tree
final.df<-d[match(mytree.names, d$taxa),]
namelist2<-final.df$taxa
namelist2==mytree.names # check 
final.df$taxa== mytree.names #check
setdiff(mytree.names,namelist2)

### try pic but im not sure how thisd work
fls.pic<-pic(final.df$meanFLS,pruned.by.anthy,,var.contrasts = TRUE,scaled = TRUE)
lm(fls.pic~1)

instread just calculate a phyosignal
rownames(final.df)<-final.df$taxa
phylosig(pruned.by.anthy,final.df$meanFLS,method="lambda",nsim = 100) ### 0.17
