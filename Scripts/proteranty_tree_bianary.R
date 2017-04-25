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
final.df$pro[final.df$Phen.sequence == "pro"] <- 1
final.df$pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$pro[final.df$Phen.sequence== "syn"] <- 0
final.df$pro[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$pro[final.df$Phen.sequence== "ser"] <- 0 
final.df$pro[final.df$Phen.sequence== "hyst"] <- 0

###now make pollination syndrom discrete
final.df["pol"]<-NA
final.df$pol[final.df$Pollination == "wind"] <- 9
final.df$pol[final.df$Pollination == "ambo"] <- 8
final.df$pol[final.df$Pollination == "insect"] <- 7

###make factors for bianary traits
final.df$pro<-as.integer(final.df$pro)
final.df$pol<-as.integer(final.df$pol)

##subset down to relevant columns
final.df<-select_(final.df, "name", "pro", "pol")

##make comprative data object
which(final.df$name%in%pruned.by.anthy$tip.label)
which(pruned.by.anthy$tip.label%in%final.df$name)
pruned.by.anthy$node.label<-""
d<-comparative.data(pruned.by.anthy,final.df,name, na.omit = FALSE)
d
##calculate phylo.d ### this syntax makes an error
phylo.d(data = final.df,phy = pruned.by.anthy, names.col = name, binvar = pro, permut = 1000, rnd.bias = NULL)
##try another way, comparative data objects should be able to just run

phylo.d(d, binvar=pro)

###try with phylobase http://www2.hawaii.edu/~mbutler/PDFs/Ch10.Phylobase.pdf
###only works as continuous with column labels
library(tidyverse)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
head(final.df)
#####picante
#https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf
final.df$pro<-as.integer(final.df$pro)
final.df$pol<-as.integer(final.df$pol)
par(mfrow=c(1,2))
 for (i in names(final.df)) {
   plot(pruned.by.anthy, show.tip.label=FALSE, main=c(i))
   tiplabels(pch=22, col=final.df[,i]+1, bg=final.df[,i]+1, cex=0.5)
 }
###signal for traits with picante
final.df <- final.df[pruned.by.anthy$tip.label,]
View(multiPhylosignal(final.df,pruned.by.anthy))



