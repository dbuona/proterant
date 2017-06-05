
#search_treebase('Acer', by='taxon')
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
#load packages
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(phylobase)
setwd("~/Documents/git/proterant/Input")
##read in tree, and data
maple<-read.nexus("T5304.nex")
names.intree<-maple$tip.label


map.dat<-read.csv("acer_data.csv")
map.dat$name<-paste("Acer",map.dat$specie_source,sep="_")
namelist<-unique(map.dat$name)

###prune tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(maple,to.prune)
mytree.names<-pruned.by.anthy$tip.label
setdiff(namelist,mytree.names)

#order them
final.df<-map.dat[match(mytree.names, map.dat$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names

library(tidyverse)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
head(final.df)
pro.df<-select_(final.df, "anthy")
colors2<-setNames(c("red","orange","yellow" ,"light blue","green","dark green"),c("hysteranthous","synanthous/hysteranthous","hysteranthous/synanthous","synanthous","synanthous/seranthous", "seranthous"))
dotTree(pruned.by.anthy,pro.df,legend=TRUE,fsize=.8,ftype="i",colors=colors2)




