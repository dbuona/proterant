##basic exploration of harvard forest data
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(plotrix)
library(gdata)
library(nlme)
library(scales)
library(arm)
library(picante)
library(ade4)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/git/hysterant/Data")
hf<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
head(hf)
reduced<-select(hf,species,l75.jd, fopn.jd)
grouped<-group_by(reduced,species)
head(grouped)

View(grouped)
sum(is.na(grouped$fopn.jd))
#how many are not NA in each catory
aggregate(fopn.jd ~ species, data=grouped, function(x) {sum(!is.na(x))}, na.action = NULL)
ag<-aggregate(l75.jd ~ species, data=grouped, function(x) {sum(!is.na(x))}, na.action = NULL)
ag$l75.jd
