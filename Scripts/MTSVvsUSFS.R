
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")

###libraryies
library("tidyverse")
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

together<-intersect(mich.data$name,silv.data$name)
mich.match<-filter(mich.data, name %in% together)
fs.match<-filter(silv.data, name %in% together)
mich.match<-dplyr::select(mich.match,name,Phen.sequence,pro,pro3,pro2)
fs.match<-dplyr::select(fs.match,name,silvic_phen_seq,pro,pro3,pro2)
goo<-data.frame(fs.match$pro3,mich.match$pro3)
goo$matcher<-ifelse(goo$fs.match.pro3==goo$mich.match.pro3,1,0)
table(goo$matcher) #30% are different emprically


### For Final plot
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
}
extract_CIs<-function(x){
  filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
}  
###read in data
mich.data<-read.csv("datasheets_derived/MTSV_USFS/michdata_final.csv")
mich.tre<-read.tree("datasheets_derived/MTSV_USFS/michtre_final.tre")

silv.data<-read.csv("datasheets_derived/MTSV_USFS/silvdata_final.csv")
silv.tre<-read.tree("datasheets_derived/MTSV_USFS/silvtre_final.tre")
##### phylo signal for suppliment
mich.tre$node.label<-NULL
silv.tre$node.label<-NULL
