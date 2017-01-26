##basic exploration of harvard forest data
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(plotrix)
library(gdata)
library(nlme)
library(scales)
library(arm)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/git/hysterant/Data")
hf<-read.csv("hf003-05-mean-ind.csv",header=TRUE)

###plot the budburst data###
new<-rename(hf,lbb.jd=bb.jd)
new<-select(new,-l75.jd)
head(new)
new<-gather(new,phenophase,eventday,lbb.jd:fbb.jd)
###Filter out species with insufficient flowering records
new<-filter(new, species %in% c( "ACPE","ACRU", "AMSP","BEAL","BEPA",
"BEPO","COAL","FRAM","ILVE","KAAN","KALA","LYLI","NEMU",
"POTR","QURU","QUVE","RHSP","SAPU","VACO","VIAL"))
q<-ggplot(new, aes(x=year, y=eventday, color=phenophase)) +
  stat_summary()+labs(title="Flower and Leaf Budburst at Harvard Forest", x="Year", y="Days since initiation")
q+facet_wrap(~species)
##Now plot flowers open vs. l75
new2<-select(hf,-fbb.jd)
head(new2)
new2<-gather(new2,phenophase,eventday,l75.jd:fopn.jd)
###Filter out species with insufficient flowering records
new2<-filter(new2, species %in% c( "ACPE","ACRU", "AMSP","BEAL","BEPA",
                                  "BEPO","COAL","FRAM","ILVE","KAAN","KALA","LYLI","NEMU",
                                  "POTR","QURU","QUVE","RHSP","SAPU","VACO","VIAL"))
q<-ggplot(new2, aes(x=year, y=eventday, color=phenophase)) +
  stat_summary()+labs(title="Flowering and Leaf Expansion at Harvard Forest", x="Year", y="Days since initiation")
q+facet_wrap(~species)

###Now, how many fopn, before lbb.
new3<-rename(hf,lbb.jd=bb.jd)
head(new3)
new3<-select(new3,-l75.jd)
new3<-select(new3,-fbb.jd)
new3<-gather(new3,phenophase,eventday,lbb.jd:fopn.jd)
new3<-filter(new3, species %in% c( "ACPE","ACRU", "AMSP","BEAL","BEPA",
                                   "BEPO","COAL","FRAM","ILVE","KAAN","KALA","LYLI","NEMU",
                                   "POTR","QURU","QUVE","RHSP","SAPU","VACO","VIAL"))
q<-ggplot(new3, aes(x=year, y=eventday, color=phenophase)) +
  stat_summary()+labs(title="Open Flowers and Leaf Budbust", x="Year", y="Days since initiation")
q+facet_wrap(~species)
