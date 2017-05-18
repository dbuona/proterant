###first attempt at plgs
## phylogeny of proterany as a bianary trait, adapt from first phylo_script.R april 10
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
#load packages
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(phylobase)
library(tidyverse)
setwd("~/Documents/git/proterant/input")

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
#list of species in tree
names.intree<-treee$tip.label

### read in data from michigan trees, format it like Zanne
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy$name<-paste(anthy$Genus,anthy$Species,sep="_")
# list of my species myspecies
namelist<-unique(anthy$name)
##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

#put the tree and data sheet in same order
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names
##okay all catagories match
#now make a new column for proteranthy as a bianary trait
final.df["pro"]<-NA
final.df$pro[final.df$Phen.sequence == "pro"] <- 1
final.df$pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$pro[final.df$Phen.sequence== "syn"] <- 0
final.df$pro[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$pro[final.df$Phen.sequence== "ser"] <- 0 
final.df$pro[final.df$Phen.sequence== "hyst"] <- 0

###now make pollination syndrom discrete
final.df["pol"]<-0
final.df$pol[final.df$Pollination == "wind"] <- 9
final.df$pol[final.df$Pollination == "ambo"] <- 8
final.df$pol[final.df$Pollination == "insect"] <- 7
final.df$pol[final.df$Pollination ==0 ] <-8

final.df$pol<-ifelse(final.df$pol==0, 8, final.df$pol)

###make factors for bianary traits
final.df$pro<-as.integer(final.df$pro)
final.df$pol<-as.integer(final.df$pol)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
head(final.df)
name.check(pruned.by.anthy, final.df)

#model
pro <- final.df[, "pro"]
pol<- final.df[, "pol"]
height<-final.df[,"heigh_height"]

pglsModel2 <- gls(pro ~ pol+height, correlation = corBrownian(phy = pruned.by.anthy),
                  data = final.df, method = "ML")
anova(pglsModel2)
###I think it needs to be a logistical regression as bianary

library("phylolm")
mod3<-phyloglm(pro~pol,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot = 0, full.matrix = TRUE)
summary(mod3)
vcov(mod3)
coef(mod3)
