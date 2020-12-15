rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/investment")

library(stringr)
library("ncdf4")
library(raster)
library(ggplot2)
library("brms")




palmer.b <- brick("lbda-v2_kddm_pmdi_2017.nc")
sp<-read.csv("Input/midwest_round1Dec11.csv") ## all the prunus records

lonpoints<-sp$lon # make vector of prunus coordinates
latpoints<-sp$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
ext<-raster::extract(mean.prunus,extract.pts,method="simple")

sp$pdsi<-ext

ggplot(sp,aes(specificEpithet,pdsi))+geom_boxplot()
ggpubr::ggboxplot(sp,'specificEpithet','pdsi')

coef(lme4::lmer(pdsi~bbch.v+(bbch.v|specificEpithet),data=sp))
           