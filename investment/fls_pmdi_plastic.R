rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
library(raster)
library(ncdf4)
library(rgdal)
library(tidyverse)
library(brms)

setwd("~/Documents/git/proterant/investment/Input")
d.flo<-read.csv("input_clean/FLS_clean.csv")
d.um<-d.flo#<-filter(d.flo,specificEpithet %in% c("maritima","mexicana","alleghaniensis","nigra","umbellata"))
palmer.b <- brick("..//Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d.um$lon # make vector of prunus coordinates
latpoints<-d.um$lat #
extract.pts <- cbind(lonpoints,latpoints)
palmer.b
palmer.b2 <-palmer.b[[1900:2018]]## subset to pnly last century
palmer.b2<-brick(palmer.b2)
ext<-raster::extract(palmer.b2,extract.pts,method="simple")
colnames(ext)<-(c(1899:2017))
ext<-as.data.frame(ext)
ext$lat<-latpoints
ext$lon<-lonpoints
goo<-tidyr::gather(ext,"year","pdsi",1:119)
class(goo$year)
g
goo$year<-as.integer(goo$year)

d.um<-left_join(d.um,goo)
ggplot(d.um,aes(pdsi,bbch.v.scale))+geom_point()+geom_smooth(method="lm")+facet_wrap(~specificEpithet)
ggplot(d.um,aes(pdsi,bbch.short))+geom_point()+geom_smooth(method="lm")+facet_wrap(~specificEpithet)
moda<-brm(bbch.v.scale~pdsi+doy.cent+(pdsi+doy.cent|specificEpithet),data=d.um,family=cumulative("logit"))
pdsiflo<-brm(doy~pdsi+lat+(pdsi+lat|specificEpithet),data=d.um)
summary(moda)
coef(moda,probs = c(.25,.75))
fixef(moda,probs = c(.25,.75))

coef(pdsiflo,probs = c(.25,.75))
fixef(pdsiflo,probs = c(.25,.75))
summary(pdsiflo)

## take away, for the most part once you control for day of the ear drier conditions seem to trigger stronger hysteranthy
## despite later flower in general do to dryness