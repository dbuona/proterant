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

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
##did it work? Let's check
pruned.by.anthy
plot(pruned.by.anthy)
mytree.names<-pruned.by.anthy$tip.label
mytree.names<-pruned.by.anthy$tip.label
setdiff(namelist,mytree.names) 
### 108 of my species were in the tree. 
#Add in the remaining ones, need to figure out how to do this with out dropping banch lengths, perhaps with pez function congeneric.merge
#species<-c("Acer_barbatum","Acer_nigrum","Aesculus_octandra","Carya_aquatica","Carya_illinoiensis","Carya_laciniosa","Celtis_tenuifolia","Fraxinus_pensylvanica","Gymnocladus_dioicus","Malus_coronaria","Populus_heterophylla","Prunus_nigra","Quercus_prinoides","Quercus_coccinea","Quercus_ellipsoidalis","Quercus_douglasii","Quercus_muehlenbergii","Quercus_nuttallii","Quercus_phellos","Quercus_prinus","Salix_nigra","Sorbus_americana","Tilia_heterophylla","Ulmus_thomasii")
#  for(i in 1:length(species)) pruned.by.anthy<-add.species.to.genus(pruned.by.anthy,species[i],where="root")
# plotTree(pruned.by.anthy,ftype="i")

 ##this step is also dealing with removing or adding species not on zanne tree. deal with thisa gain in the future
 #for now, I will subset to exclude them
#anthy<-filter(anthy,name != "Gymnocladus_dioicus")
 #anthy<-filter(anthy,name !="Populus_heterophylla")
 #anthy<-filter(anthy,name !="Prunus_nigra") 
 #anthy<-filter(anthy,name !="Quercus_prinoides")
 #anthy<-filter(anthy,name !="Quercus_coccinea")
 #anthy<-filter(anthy,name !="Quercus_ellipsoidalis" )
 #anthy<-filter(anthy,name != "Quercus_douglasii")
 #anthy<-filter(anthy,name != "Quercus_muehlenbergii")
 #anthy<-filter(anthy,name != "Quercus_nuttallii")
 #anthy<-filter(anthy,name !="Quercus_phellos")
 #anthy<-filter(anthy,name !="Quercus_prinus")
 #anthy<-filter(anthy,name != "Salix_nigra")
 #anthy<-filter(anthy,name !="Sorbus_americana")
 #anthy<-filter(anthy,name !="Tilia_heterophylla")
# anthy<-filter(anthy,name !="Ulmus_thomasii")
##species list
 
namelist<-unique(anthy$name)
setdiff(namelist,mytree.names) 
###okay, now there is no differentce between namelist and mytree.names

final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names
###matches!
final.df["pro"]<-NA

  
final.df$pro[final.df$Proteranthous == 1] <- 3
final.df$pro[final.df$Proteranthous == 0] <- 2
final.df$pro[final.df$Proteranthous== NA] <- 1
  
#which(anthy$name%in%pruned.by.anthy$tip.label)
#which(pruned.by.anthy$tip.label%in%anthy$name)
#pruned.by.anthy$node.label<-""
#comparative.data(pruned.by.anthy,anthy,names.col=name,na.omit=FALSE)

 ##plot on tree
plot.phylo(pruned.by.anthy,show.tip.label = TRUE, tip.color= final.df$pro, adj=.4,align.tip.label = TRUE, cex = 0.5)
###this worked! green=proteranthous, red= non-proteranthous, black= NA
###method below doesnt work
x<-final.df$Proteranthous
names(x)=pruned.by.anthy$tip.label

phylosig(pruned.by.anthy, x, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())
help("phylosig")

