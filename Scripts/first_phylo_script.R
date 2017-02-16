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
library(picante)

setwd("~/Documents/git/proterant/input")

###read in tree
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
names.intree<-treee$tip.label

### read in data, format it like Zanne
anthy<-read.csv("proterant_ds.csv", header = TRUE)
anthy$name<-paste(anthy$genus,anthy$species,sep="_")

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
 which(anthy$name%in%pruned.by.anthy$tip.label)
 which(pruned.by.anthy$tip.label%in%anthy$name)
 pruned.by.anthy$node.label<-""
 signal<-comparative.data(pruned.by.anthy,anthy,names.col=name,na.omit=FALSE)
signal


##############This was troubleshooting an old problem, but it works now###############
#anthy$name
#troubleshoot
#intersect(pruned.by.anthy$node.label,pruned.by.anthy$tip.label)
#name.check(phy = pruned.by.anthy,data = anthy,data.names = anthy$name)
#pruned.by.anthy$node.label<-""
#pruned.by.anthy$tip.label
######################################################################################
#follow ucdavis workflow for discrete traits
###new column for bianry proteranthy or non proteranthy
anthy$proteranthy[anthy$mich_phen_seq == "pro"] <- 1
anthy$proteranthy[anthy$mich_phen_seq == "pro/syn"] <- 1
anthy$proteranthy[anthy$mich_phen_seq == "syn"] <- 0
anthy$proteranthy[anthy$mich_phen_seq == "syn/ser"] <- 0
anthy$proteranthy[anthy$mich_phen_seq == "ser"] <- 0
anthy$proteranthy[anthy$mich_phen_seq == "hyst"] <- 0
View(anthy)

### trying to make a plot where tip labels are color coded if proteranthous or not, non below work

plot.phylo(pruned.by.anthy,show.tip.label = TRUE, tip.color = "blue")
trait.plot(pruned.by.anthy,dat = anthy, cols = anthy$proteranthy)   
       
###markov####these won't run
help(fitDiscrete)
pro_equalrate<-fitDiscrete(pruned.by.anthy$proteranty,dat=anthy, model="ER")
pro_ardrate<-fitDiscrete(pruned.by.anthy,dat  anthy$proteranthy, model="ARD")
###lamda


