###THis compares phyloD with differnt tree addition methods
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(tidyverse)
library(boot)
library(phylolm)
#https://academic-oup-com.ezp-prod1.hul.harvard.edu/sysbio/article-lookup/doi/10.1093/sysbio/syp074 Garland and Ives 2010

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy<-filter(anthy,Phen.sequence!="evergreen")
anthy<-filter(anthy,Phen.sequence!="non_woody")
anthy<-filter(anthy,Phen.sequence!="unknown")

names.intree<-treee$tip.label

#dataformat it like Zanne
anthy$name<-paste(anthy$Genus,anthy$Species,sep="_")

# list of my species myspecies
namelist<-unique(anthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
#plot(pruned.by.anthy)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #158 species include

addins<-setdiff(namelist,mytree.names)

is.ultrametric(pruned.by.anthy)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)

###phyloD with no additions. Based on 158 species
###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names


setdiff(mytree.names,namelist2)


####add comlumns for analysis
final.df["pro"]<-NA
final.df$pro[final.df$Phen.sequence == "pro"] <- 1
final.df$pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$pro[final.df$Phen.sequence== "syn"] <- 0
final.df$pro[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$pro[final.df$Phen.sequence== "ser"] <- 0 
final.df$pro[final.df$Phen.sequence== "hyst"] <- 0

pruned.by.anthy$node.label<-NULL
final.df<-na.omit(final.df)
plot(pruned.by.anthy)

zoom(pruned.by.anthy, grep("Salix", pruned.by.anthy$tip.label))
zoom(pruned.by.anthy, grep("Acer", pruned.by.anthy$tip.label))
zoom(pruned.by.anthy, grep("Rhus", pruned.by.anthy$tip.label))
zoom(pruned.by.anthy, grep("Prunus", pruned.by.anthy$tip.label))
zoom(pruned.by.anthy, grep("Vaccinium", pruned.by.anthy$tip.label))
###Acer, Rhus, Salix, Prunus,Vaccinium, Sassafras.Lindera
#Rutaceae: Znthoxylem/Ptelea
#Corylus. Ostrya/Caprinus, 
#Eleagnaceae
##
##Other possibilities from different trees



