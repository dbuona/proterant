rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/")

library(ape)
library(phytools)
library(brms)
library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")

HF<-read.csv("Data/hf003-05-mean-ind.csv",header=TRUE)

### This file is for prepping USDA traits for merging with hf
traits<-read.csv("Input/HarvardForest/HF.species.trait.data.csv",header=TRUE)
traits<-dplyr::filter(traits,species!="QUAL")
###make things numeric
traits$pol<-ifelse(traits$pollination=="insect",0,1)
traits$flo_view<-ifelse(traits$fl_conspic=="Y",0,1)
traits$frost_free<-as.numeric(traits$frost_free)

###prep the tree
treee<-read.tree("Data/Vascular_Plants_rooted.dated.tre")
names.intree<-treee$tip.label
namelist<-unique(traits$name)

to.prune<-which(!names.intree%in%namelist)
HF.tree.pruned<-drop.tip(treee,to.prune)
mytree.names<-HF.tree.pruned$tip.label

setdiff(namelist,mytree.names) ###Viburum
intersect(namelist,mytree.names)

###make ultrametric
HF.tree.pruned<-chronoMPL(HF.tree.pruned)
is.ultrametric(HF.tree.pruned)

###add one species
HF.tree.pruned<-add.species.to.genus(HF.tree.pruned, "Viburnum_lantanoides",genus=NULL,where="root")

##matchthem ### not sure if I need to do this for brms
df<-traits[match(mytree.names, traits$name),]
namelist<-df$name
namelist==mytree.names
df$name== mytree.names
#### cool
HF.tree.pruned$node.label<-NULL

write.csv(df,"sub_projs/HarvardForest/HFdata4modeling.csv")
write.tree(HF.tree.pruned,"sub_projs/HarvardForest/HFtree4modeling.tre")
