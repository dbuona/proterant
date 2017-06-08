###This will be the source code for calling my phylogenetic tree adn michigan trees data to be used in further analysis
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
#load packages

library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)

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
#plot(pruned.by.anthy)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #80 species include

setdiff(namelist,mytree.names) #14 species did not make it

##Add missing species
#[1] "Populus_heterophylla"  "Salix_nigra"           "Gymnocladus_dioicus"   "Malus_coronaria"       "Sorbus_americana"      "Prunus_nigra"         
#[7] "Ulmus_thomasii"        "Celtis_tenuifolia"     "Quercus_muehlenbergii" "Quercus_prinoides"     "Quercus_coccinea"      "Quercus_ellipsoidalis"
#[13] "Carya_lacinosa"        "Acer_nigrum" 
###or not

###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names

stop("this is not an error, just tells you we've finished the code")





