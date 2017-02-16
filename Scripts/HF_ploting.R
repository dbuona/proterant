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

setwd("~/Documents/git/proterant/Data")
hf<-read.csv("hf003-05-mean-ind.csv",header=TRUE)

###plot the budburst data###
new<-rename(hf,lbb.jd=bb.jd)
head(new)
new<-dplyr::select(new,-l75.jd)
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
new2<-dplyr::select(hf,-fbb.jd)
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
new3<-dplyr::select(new3,-l75.jd)
new3<-dplyr::select(new3,-fbb.jd)
new3<-gather(new3,phenophase,eventday,lbb.jd:fopn.jd)
new3<-filter(new3, species %in% c( "ACPE","ACRU", "AMSP","BEAL","BEPA",
                                   "BEPO","COAL","FRAM","ILVE","KAAN","KALA","LYLI","NEMU",
                                   "POTR","QURU","QUVE","RHSP","SAPU","VACO","VIAL"))
q<-ggplot(new3, aes(x=year, y=eventday, color=phenophase)) +
  stat_summary()+labs(title="Open Flowers and Leaf Budbust", x="Year", y="Days since initiation")
q+facet_wrap(~species)
View(hf)

### now try to make an anual mean leaf out value
mean<-read.csv("hf003-06-mean-spp.csv",header=TRUE)
head(mean)
mean<-rename(mean,lbb.jd=bb.jd)
head(mean)
#reduced_mean<-gather(mean,phenophase,eventday,lbb.jd:fopn.jd)
reduced_mean<-filter(mean, species %in% c( "ACPE","ACRU", "AMSP","BEAL","BEPA",
                                 "BEPO","COAL","FRAM","ILVE","KAAN","KALA","LYLI","NEMU",
                                 "POTR","QURU","QUVE","RHSP","SAPU","VACO","VIAL"))
head(reduced_mean)
reduced_mean<-group_by(reduced_mean, year)
leafout<-summarise(reduced_mean, avgleaf=mean(l75.jd,na.rm=TRUE), sdleaf=sd(l75.jd, na.rm=TRUE))
leafout
flo_point<-group_by(reduced_mean,year, species)
floers<-summarise(flo_point, avgfl=mean(fopn.jd,na.rm=TRUE), sdfl=sd(fopn.jd, na.rm=TRUE))
floers
##Average annual flower open date per species compared to average 75 leaf size in community
r<-ggplot(floers, aes(x=year, y=avgfl)) +geom_point(aes(x=year, y=avgfl,color=species))+geom_smooth(data = leafout, aes(x=year, y=avgleaf, colour="Average Leafout"))
r


