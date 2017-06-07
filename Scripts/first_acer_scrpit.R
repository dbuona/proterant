
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
#mapl<-read.nexus("T5305.nex") ## no branch lengths
maple<-read.nexus("T5304.nex") ##https://treebase.org/treebase-web/search/treeSearch.html oe "T5305.nex" from Renner et al 2008 no branch lengths
maple$edge.length
#maple<- #how to choose what branch lengths to assign. See below.
#This choice definitely effects results
plot(compute.brlen(maple, runif, min = 0, max = 5))
layout(matrix(1:4, 2, 2))
plot(compute.brlen(maple, power=1), main=expression(rho==1))
plot(compute.brlen(maple, power=3), main=expression(rho==3))
plot(compute.brlen(maple, power=0.5), main=expression(rho==0.5))
plot(compute.brlen(maple, power=0.1), main=expression(rho==0.1))
layout(1)
maple<-compute.brlen(maple, power=0.5)

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
colors2<-setNames(c("red","red","red" ,"light blue","dark green","dark green"),c("hysteranthous","synanthous/hysteranthous","hysteranthous/synanthous","synanthous","synanthous/seranthous", "seranthous"))
dotTree(pruned.by.anthy,pro.df,legend=TRUE,fsize=.8,ftype="i",colors=colors2)

###signal
library("phylolm")
final.df$seq<-0
final.df$seq[final.df$name == "synanthous"] <- 0
final.df$seq[final.df$name == "synanthous/seranthous"] <- 0
final.df$seq[final.df$name == "seranthous"] <- 0
final.df$seq[final.df$anthy == "hysteranthous"] <- 1
final.df$seq[final.df$anthy == "hysteranthous/synanthous"] <- 1
final.df$seq[final.df$name == "synanthous/hysteranthous"] <- 1

###physignal
final.df$seq<-as.integer(final.df$seq)
pro.df<-select_(final.df, "seq")
for (i in names(pro.df)) {
  plot(pruned.by.anthy, show.tip.label=FALSE, main=c(i))
  tiplabels(pch=22, col=final.df[,i]+1, bg=final.df[,i]+1, cex=0.5)
}
pro.df$seq<-as.character(pro.df$seq)
colors3<-setNames(c("pink", "dark green"),c("1","0"))
dotTree(pruned.by.anthy,pro.df,legend=TRUE,fsize=.8,ftype="i",colors=colors3)
##phylogenetic signal
phylosig(pruned.by.anthy,final.df$seq,test=TRUE)
phylosig(pruned.by.anthy,final.df$seq,method="lambda",test=TRUE)

##future model
mod4a<-phyloglm(seq~1,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod4a)
#library(nnet)
#multinom(anthy~1,data=final.df)


