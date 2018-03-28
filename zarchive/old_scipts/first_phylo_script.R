##DAn B start phylogenies on 1/19/2017
#adapted for michigan trees data only on 3/12
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
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy$name<-paste(anthy$Genus,anthy$Species,sep="_")

###myspecies
namelist<-unique(anthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
##did it work? Let's check
pruned.by.anthy
plot(pruned.by.anthy)
mytree.names<-pruned.by.anthy$tip.label

setdiff(namelist,mytree.names) 
intersect(namelist,mytree.names)
### 80 of my species were in the tree. 
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

####now for functional
final.df["func.pro"]<-NA
final.df$func.pro[final.df$Phen.sequence == "pro"] <- 1
final.df$func.pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$func.pro[final.df$Phen.sequence== "syn"] <- 2
final.df$func.pro[final.df$Phen.sequence== "syn/ser"] <- 2
final.df$func.pro[final.df$Phen.sequence== "ser"] <- 3  
final.df$func.pro[final.df$Phen.sequence== "hyst"] <- 3 
##proteranthous as bianary
final.df[".pro"]<-NA
final.df$bin.pro[final.df$Phen.sequence == "pro"] <- 1
final.df$bin.pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$bin.pro[final.df$Phen.sequence== "syn"] <- 0
final.df$bin.pro[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$bin.pro[final.df$Phen.sequence== "ser"] <- 0  
final.df$bin.pro[final.df$Phen.sequence== "hyst"] <- 0 

 ##plot on tree
plot.phylo(pruned.by.anthy,show.tip.label = TRUE, tip.color= final.df$func.pro, adj=.4,align.tip.label = TRUE, cex = 0.5)
plot.phylo(pruned.by.anthy,show.tip.label = TRUE, tip.color= final.df$pro, adj=.4,align.tip.label = TRUE, cex = 0.5)
plot.phylo(pruned.by.anthy,show.tip.label = TRUE, tip.color= final.df$bin.pro, adj=.4,align.tip.label = TRUE, cex = 0.5)
###this worked! green=proteranthous, red= non-proteranthous, black= NA

##try plotting with picante
traits<-dplyr::select(final.df,name,pro,bin.pro)
#par(mfrow=c(2,2))
 for (i in names(traits)) {
   plot(anthy, show.tip.label=FALSE, main=i)
   tiplabels(pch=22, col=traits[,i]+1, bg=traits[,i]+1, cex=1.5) }


###preparing to calculate lamda
#x<-final.df$pro #something seems wrong
#names(x)=pruned.by.anthy$tip.label

#y<-final.df$func.pro
#names(y)=pruned.by.anthy$tip.label

#z<-final.df$bin.pro
#names(z)=pruned.by.anthy$tip.label
####lamda test- these are all the same
#phylosig(pruned.by.anthy, x, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
#         control=list())
#phylosig(pruned.by.anthy, y, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
#         control=list())
#phylosig(pruned.by.anthy, z, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
#         control=list())
###nacho says not to use lamda because its for continuous traits
which(final.df$name%in%pruned.by.anthy$tip.label)
which(pruned.by.anthy$tip.label%in%final.df$name)
pruned.by.anthy$node.label<-""

d<-comparative.data(pruned.by.anthy,final.df,name, na.omit = FALSE)
phylo.d(d$final.df,pruned.by.anthy,names.col = name ,binvar =bin.pro ,permut = 1000, rnd.bias=NULL)

