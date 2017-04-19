## phylogeny of proterany as a bianary trait, adapt from first phylo_script.R april 10
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
#load packages
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(randomForest)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(phylobase)
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
final.df$pro[final.df$Phen.sequence == "pro"] <- "pro"
final.df$pro[final.df$Phen.sequence == "pro/syn"] <- "pro"
final.df$pro[final.df$Phen.sequence== "syn"] <- "syn"
final.df$pro[final.df$Phen.sequence== "syn/ser"] <- "ser"
final.df$pro[final.df$Phen.sequence== "ser"] <- "ser" 
final.df$pro[final.df$Phen.sequence== "hyst"] <- "hyst"

###now make pollination syndrom discrete
final.df["pol"]<-NA
final.df$pol[final.df$Pollination == "wind"] <- "wind"
final.df$pol[final.df$Pollination == "ambo"] <- "ambo"
final.df$pol[final.df$Pollination == "insect"] <- "insect"
final.df<-select_(final.df, "name", "pro", "pol")
library(tidyverse)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
head(final.df)
final.df$pro<-as.factor(final.df$pro)
final.df$pol<-as.character(final.df$pol)
colors<-setNames(c("pink","yellow","green","grey","light blue", "navy blue","orange" ),c("pro","syn","ser","hyst","wind","ambo","insect"))
dotTree(pruned.by.anthy,final.df,colors=colors,data.type="discrete",fsize=0.7,x.space=0.05)

###plot for hysteranthy only
pro.df<-select_(final.df, "pro")
colors2<-setNames(c("red","yellow","dark green","orange" ),c("pro","syn","ser","hyst"))
dotTree(pruned.by.anthy,pro.df,legend=TRUE,fsize=.5,ftype="i",colors=colors2)

