
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")

###libraryies
library(dplyr)
library("ggplot2")

library(brms)
library(ggthemes)
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library(boot)
library("phylolm")
library(ggstance)

library("raster")
library("remote")
library(reshape2)
library(RColorBrewer)
mich.data<-read.csv("..//sub_projs/MTSV_USFS/mich_data_full_clean.csv")
silv.data<-read.csv("..//sub_projs/MTSV_USFS/silv_data_full.csv")
together<-intersect(mich.data$name,silv.data$name)
mich.match<-dplyr::filter(mich.data, name %in% together)
fs.match<-dplyr::filter(silv.data, name %in% together)
mich.match<-dplyr::select(mich.match,name,Phen.sequence)
fs.match<-dplyr::select(fs.match,name,silvic_phen_seq)
goo<-data.frame(fs.match$silvic_phen_seq,mich.match$Phen.sequence)
goo$matcher<-ifelse(goo$fs.match.silvic_phen_seq==goo$mich.match.Phen.sequence,1,0)
table(goo$matcher) #30% are different emprically
15/(49) ### 30%
write.csv(mich.data,"../sub_projs/MTSV_USFS/michdata_final.csv")
write.tree(mich.tree.droughtprune,"../sub_projs//MTSV_USFS/michtre_final.tre")

