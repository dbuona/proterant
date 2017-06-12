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


###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
fullanthy<-read.csv("proterant_ds.csv", header = TRUE)

names.intree<-treee$tip.label

#dataformat it like Zanne
fullanthy$name<-paste(fullanthy$genus,fullanthy$species,sep="_")

# list of my species myspecies
namelist<-unique(fullanthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.fullanthy<-drop.tip(treee,to.prune)
#plot(pruned.by.anthy)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.fullanthy$tip.label

intersect(namelist,mytree.names) 
setdiff(namelist,mytree.names)

final.df2<-fullanthy[match(mytree.names, fullanthy$name),]
namelist2<-final.df2$name
namelist2==mytree.names
final.df2$name== mytree.names

final.df2["michpro"]<-NA
final.df2["silvpro"]<-NA

final.df2$michpro[final.df2$mich_phen_seq == "pro"] <- 1
final.df2$michpro[final.df2$mich_phen_seq == "pro/syn"] <- 1
final.df2$michpro[final.df2$mich_phen_seq== "syn"] <- 0
final.df2$michpro[final.df2$mich_phen_seq== "syn/ser"] <- 0
final.df2$michpro[final.df2$mich_phen_seq== "ser"] <- 0 
final.df2$michpro[final.df2$mich_phen_seq== "hyst"] <- 0

final.df2$silvpro[final.df2$silvic_phen_seq == "pro"] <- 1
final.df2$silvpro[final.df2$silvic_phen_seq == "pro/syn"] <- 1
final.df2$silvpro[final.df2$silvic_phen_seq== "syn"] <- 0
final.df2$silvpro[final.df2$silvic_phen_seq== "syn/ser"] <- 0
final.df2$silvpro[final.df2$silvic_phen_seq== "ser"] <- 0 
final.df2$silvpro[final.df2$silvic_phen_seq== "hyst"] <- 0

#### make intergrated response which favors michigan trees
final.df2["proM"]<-NA
final.df2$proM <-ifelse(final.df2$michpro==final.df2$silvpro, final.df2$michpro, final.df2$proM)
final.df2$proM <-ifelse(is.na(final.df2$michpro), final.df2$silvpro, final.df2$proM) 
final.df2$proM <-ifelse(is.na(final.df2$silvpro), final.df2$michpro, final.df2$proM)

final.df2$shrink<-NA 
final.df2$shrink<-ifelse(is.na(final.df2$michpro)& is.na(final.df2$silvpro), 1,0)
final.df2<-filter(final.df2, shrink==0)

final.df2$dum<-NA
final.df2$proM<-ifelse(is.na(final.df2$proM) , final.df2$michpro, final.df2$proM)
##make integrated resposne column that favors silvics

final.df2$proS<-NA
final.df2$proS <-ifelse(final.df2$michpro==final.df2$silvpro, final.df2$michpro, final.df2$proS)
final.df2$proS <-ifelse(is.na(final.df2$michpro), final.df2$silvpro, final.df2$proS) 
final.df2$proS <-ifelse(is.na(final.df2$silvpro), final.df2$michpro, final.df2$proS)

final.df2$shrink<-NA 
final.df2$shrink<-ifelse(is.na(final.df2$michpro)& is.na(final.df2$silvpro), 1,0)
final.df2<-filter(final.df2, shrink==0)


final.df2$proS<-ifelse(is.na(final.df2$proS) , final.df2$silvpro, final.df2$proS)

##okay proS and proM are the response varaible of choice now.

###assign ambo
final.df2$mich_pollin[final.df2$name == "Acer_spicatum"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Acer_pensylvanicum"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Acer_negundo"] <- "wind"
final.df2$mich_pollin[final.df2$name == "Acer_platanoides"] <- "wind"
final.df2$mich_pollin[final.df2$name == "Acer_rubrum"] <- "wind"
final.df2$mich_pollin[final.df2$name == "Acer_saccharinum"] <- "wind"
final.df2$mich_pollin[final.df2$name == "Acer_saccharum"] <- "wind"
final.df2$mich_pollin[final.df2$name == "Populus_nigra"] <- "wind"
final.df2$mich_pollin[final.df2$name == "Fraxinus_latifolia"] <- "wind"
final.df2$silvics_pollin[final.df2$silvics_pollin == "ambo"] <- "insect"
final.df2$silvics_pollin[final.df2$silvics_pollin == ""] <- NA
final.df2$mich_pollin[final.df2$name == "Salix_alba"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Salix_amygdaloides"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Salix_fragilis"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Castanea_dentata"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Gleditsia_triacanthos"] <- "insect"
final.df2$mich_pollin[final.df2$name == "Nyssa_aquatica"] <- "insect"
final.df2$mich_pollin[final.df2$mich_pollin == "wnd"] <- "wind"

###now make pollination one column

final.df2$syndrome<-NA
final.df2$syndrome <-ifelse(final.df2$mich_pollin==final.df2$silvics_pollin, final.df2$mich_pollin, final.df2$syndrome)
final.df2$syndrome <-ifelse(is.na(final.df2$mich_pollin), final.df2$silvics_pollin, final.df2$syndrome) 
final.df2$syndrome <-ifelse(is.na(final.df2$silvics_pollin), final.df2$mich_pollin, final.df2$syndrome)
final.df2$shrink2<-NA 
final.df2$shrink2<-ifelse(is.na(final.df2$mich_pollin)& is.na(final.df2$silvics_pollin), 1,0)

final.df2<-filter(final.df2, shrink2==0)
