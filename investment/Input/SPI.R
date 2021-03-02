## extrac spi
### Explore the prunus data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library(tibble)
library(lubridate)
library(stringr)
library("ncdf4")
library(raster)


setwd("~/Documents/git/proterant/investment/")
spi <- brick("Data/470067.spi12.spi3_6_12_1deg_cru_ts_3_21_1949_2012.nc")
spi

d<-read.csv("input/midwest_round1Dec11.csv") # active datasheet

indices <- format(as.Date(names(spi), format = "X%Y.%m.%d"), format = "%m")
indices <- as.numeric(indices)
n<-names(spi)
nn <- as.integer(substr(n,7,8))

subj <-raster::subset(spi, which(nn %in% c(3,4,5)))
names(subj)
Monthspi<- stackApply(spi, indices, fun = mean)

minyear <-calc(spi, fun = min,na.rm=TRUE) 
minspring <-calc(subj, fun = min,na.rm=TRUE) 




lonpoints<-d$lon # make vector of prunus coordinates
latpoints<-d$lat #



extract.pts <- cbind(lonpoints,latpoints)

ext<-raster::extract(minyear,extract.pts,method="simple")
ext2<-raster::extract(minspring,extract.pts,method="simple")

d$spiyear<-ext
d$spispring<-ext2



spi.means<-d %>%group_by((specificEpithet)) %>% summarise(minspi=mean(spiyear,na.rm=TRUE),sdspi=sd(spiyear,na.rm=TRUE))
spispring.means<-d %>%group_by((specificEpithet)) %>% summarise(minspi=mean(spispring,na.rm=TRUE),sdspi=sd(spispring,na.rm=TRUE))

spi.mins<-d %>%group_by((specificEpithet)) %>% summarise(minspi=min(spiyear,na.rm=TRUE),sdspi=sd(spiyear,na.rm=TRUE))
spispring.mins<-d %>%group_by((specificEpithet)) %>% summarise(minspi=min(spispring,na.rm=TRUE),sdspi=sd(spispring,na.rm=TRUE))

quicky<-brm(spispring~(1|specificEpithet),data=d)

q<-as.data.frame(coef(quicky))
