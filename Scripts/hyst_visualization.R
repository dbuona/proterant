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
#https://academic-oup-com.ezp-prod1.hul.harvard.edu/sysbio/article-lookup/doi/10.1093/sysbio/syp074 Garland and Ives 2010

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy<-filter(anthy, !is.na(av_fruit_time))
source("source/prune_tree.R")
is.ultrametric(pruned.by.anthy)

###try tree
final.df<-dplyr::select(final.df,name,pro,pol, class2, flo_type, fruit_bin, shade_bin)

final.df<- within(final.df, pro[pro==0]<-"non-hyst")
final.df<- within(final.df, pro[pro==1]<-"hyst")
final.df<- within(final.df, pol[pol==1]<-"wind")
final.df<- within(final.df, pol[pol==0]<-"insect")
final.df<- within(final.df, class2[class2==1]<-"tree")
final.df<- within(final.df, class2[class2==0]<-"shrub")
final.df<- within(final.df, fruit_bin[fruit_bin==0]<-"early")
final.df<- within(final.df, fruit_bin[fruit_bin==1]<-"late")


library(tidyverse)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
colors<-setNames(c("pink","dark green","brown","blue","grey", "black","yellow","navy blue","light green", "orange","beige","purple" ),c("hyst","non-hyst","insect","wind","shrub","tree","unisexual", "bisexual", "early", "late","tolerant","intolerant"))

dotTree(pruned.by.anthy,final.df,data.type="discrete",colors=colors, fsize=0.7,x.space=0.05)


final.df<- within(final.df, pro[pro=="non-hyst"]<-0)
final.df<- within(final.df, pro[pro=="hyst"]<-1)
final.df<- within(final.df, pol[pol=="wind"]<-1)
final.df<- within(final.df, pol[pol=="insect"]<-0)
final.df<- within(final.df, class2[class2=="tree"]<-1)
final.df<- within(final.df, class2[class2=="shrub"]<-0)
final.df<- within(final.df, fruit_bin[fruit_bin=="early"]<-0)
final.df<- within(final.df, fruit_bin[fruit_bin=="late"]<-1)
final.df<- within(final.df, shade_bin[shade_bin=="tolerant"]<-1)
final.df<- within(final.df, shade_bin[shade_bin=="intolerant"]<-0)
final.df<- within(final.df, flo_type[flo_type=="bisexual"]<-0)
final.df<- within(final.df, flo_type[flo_type=="unisexual"]<-1)

final.df$pro<-as.numeric(final.df$pro)
final.df$pol<-as.numeric(final.df$pol)
final.df$class2<-as.numeric(final.df$class2)
final.df$flo_type<-as.numeric(final.df$flo_type)
final.df$shade_bin<-as.numeric(final.df$shade_bin)
final.df$fruit_bin<-as.numeric(final.df$fruit_bin)



par(mfrow=c(2,3))
for (i in names(final.df)) {
  plot(pruned.by.anthy, show.tip.label=TRUE, main=c(i))
  tiplabels(pch=22, col=final.df[,i]+1, bg=final.df[,i]+1, cex=0.5)
}
par(mfrow=c(1,1))
dotTree(pruned.by.anthy,final.df,data.type="discrete", fsize=0.7,x.space=0.05)
