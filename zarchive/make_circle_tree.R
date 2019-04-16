####This makes a phylogenetic tree with hysteranthy as a trait. 1/10/18
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

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("ggtree", version = "3.8")

library(ggtree)
setwd("~/Documents/git/proterant/input")


mich.data<-read.csv("datasheets_derived/MTSV_USFS/michdata_final.csv")
mich.tre<-read.tree("datasheets_derived/MTSV_USFS/michtre_final.tre")

silv.data<-read.csv("datasheets_derived/MTSV_USFS/silvdata_final.csv")
silv.tre<-read.tree("datasheets_derived/MTSV_USFS/silvtre_final.tre")


HF.data<-read.csv("HarvardForest/HF.means.data.csv")
HF.tre<-read.tree("HarvardForest/HFtree.tre")

mich.data$synth<-mich.data$pro2+mich.data$pro3
mich.data$synth2[which(mich.data$synth==0)] <- "always seranthous"
mich.data$synth2[which(mich.data$synth==1)] <- "transitional"
mich.data$synth2[which(mich.data$synth==2)] <- "always hysteranthous"

tr <-mich.tre

par(mar=c(0,0,0,0))


dd <- data.frame(taxa  = mich.data$name, hysteranthy = mich.data$synth2 )
#row.names(dd) <- NULL

jpeg("..//figure/michtreeplot.jpeg",width = 801, height = 778)
p<-ggtree(tr,layout="radial")
p <- p %<+% dd +geom_tippoint(aes(color=hysteranthy,shape=hysteranthy),size=4)
p<-p+theme(legend.position="bottom")
p
dev.off()
###run concept.paper.analysis.r
?geom_tiplab2()
silv.data$synth<-silv.data$pro2+silv.data$pro3
silv.data$synth2[which(silv.data$synth==0)] <- "always seranthous"
silv.data$synth2[which(silv.data$synth==1)] <- "transitional"
silv.data$synth2[which(silv.data$synth==2)] <- "always hysteranthous"

#### make a tree for svics

par(mar=c(.01,.01,.01,.01))
tr1 <-silv.tre

dd1 <- data.frame(taxa  = silv.data$name, hysteranthy = silv.data$synth2 )

jpeg("..//figure/silvtreeplot.jpeg",width = 801, height = 778)
q<-ggtree(tr1,layout="radial")
q <- q %<+% dd1 +geom_tippoint(aes(color=hysteranthy,shape=hysteranthy),size=4)
q<-q+theme(legend.position="bottom")
q
dev.off()

HF.data$synth<-HF.data$funct.bin+HF.data$phys.bin
HF.data$synth2[which(HF.data$synth==0)] <- "always seranthous"
HF.data$synth2[which(HF.data$synth==1)] <- "transitional"
HF.data$synth2[which(HF.data$synth==2)] <- "always hysteranthous"

dd2 <- data.frame(taxa  = HF.data$X, hysteranthy = HF.data$synth2 )
par(mar=c(.01,.01,.01,.01))
jpeg("..//figure/HFtreeplot.jpeg",width = 801, height = 778)
y<-ggtree(HF.tre,layout="radial")
y <- y %<+% dd2 +geom_tippoint(aes(color=hysteranthy,shape=hysteranthy),size=4)
y<-y+theme(legend.position="bottom")
y
dev.off()
