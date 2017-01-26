##DAn B start phylogenies on 1/19/2017
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
setwd("~/Documents/git/proterant/input")

###read in tree
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
names.intree<-treee$tip.label

### read in data, format it like Zanne
anthy<-read.csv("proterant_ds.csv", header = TRUE)
anthy$name<-paste(anthy$genus,anthy$species,sep="_")
View(anthy)
###myspecies
namelist<-unique(anthy$name)
namelist

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
##did it work? Let's check
plot(pruned.by.anthy)
mytree.names<-pruned.by.anthy$tip.label
mytree.names
### 108 of my species were in the tree. Which species didn't male it
setdiff(namelist,mytree.names)

###Add in the remaining ones
species<-c("Acer_barbatum","Acer_nigrum","Aesculus_octandra","Carya_aquatica","Carya_illinoiensis","Carya_laciniosa","Celtis_tenuifolia","Fraxinus_pensylvanica","Gymnocladus_dioicus","Malus_coronaria","Populus_heterophylla","Prunus_nigra","Quercus_prinoides","Quercus_coccinea","Quercus_ellipsoidalis","Quercus_douglasii","Quercus_muehlenbergii","Quercus_nuttallii","Quercus_phellos","Quercus_prinus","Salix_nigra","Sorbus_americana","Tilia_heterophylla","Ulmus_thomasii")
 for(i in 1:length(species)) pruned.by.anthy<-add.species.to.genus(pruned.by.anthy,species[i],where="random")
 plotTree(pruned.by.anthy,ftype="i")
 mytree.names<-pruned.by.anthy$tip.label
 setdiff(namelist,mytree.names) 
 ##forwhatever reason, 15 species could not be added (8 quercus, 1 ulmus, 1 tilia, sorbus, salix, populus and gymocladus)
 mytree.names
 #leaves me with 117 for analysis

 #seeking phylogenetic signal
#error says "Labels duplicated between tips and nodes in phylogeny"
signal<-comparative.data(pruned.by.anthy,anthy,names.col="name")

#troubleshoot
intersect(pruned.by.anthy$node.label,pruned.by.anthy$tip.label)
name.check(phy = pruned.by.anthy,data = anthy,data.names = anthy$name)
### the names match, not exactly sure what a label is and can seem to find it online


