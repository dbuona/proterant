####main cleaning file for MTSV hysteranthy dataset

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



treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
####remove taxa where hysteranthy is irrelevant or unknown.
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy<-filter(anthy,Phen.sequence!="evergreen")
anthy<-filter(anthy,Phen.sequence!="non_woody")
anthy<-filter(anthy,Phen.sequence!="unknown")

source("cleaning/prune_tree.R") 
# This file: cleans names, prunes tree to match data, adds, taxa in data to tree at root
#makes tree ultrametri
#Converts predictors to 1s and 0s
write.tree(pruned.by.anthy, "pruned_for_mich.tre")
write.csv(final.df, "mich_data_full.csv",row.names = FALSE)
source("cleaning/clean_misc.R")
write.csv(mich.data,"mich_data_full_clean.csv",row.names=FALSE)
