rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/investment")
library("housingData")
library(stringr)
library("ncdf4")
library(raster)
library(ggplot2)
library("brms")

sp<-read.csv("rhody.data.csv")

geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" ## make column names compatible

rhody.data<-dplyr::left_join(sp,geoCounty,by=c("county","stateProvince")) ## This assigns each country coordinates



palmer.b <- brick("lbda-v2_kddm_pmdi_2017.nc")  ## read in balmer drought index data 

lonpoints<-rhody.data$lon # make vector of rhody coordinates
latpoints<-rhody.data$lat #
extract.pts <- cbind(lonpoints,latpoints) #make coordinates to extract
ext <- extract(palmer.b,extract.pts,method="simple") ###extract drought info from each coordinate of prunus colllect



mean.rhody<- calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
ext <- extract(mean.rhody,extract.pts,method="simple")

rhody.data$pdsi<-ext

hyster<-read.csv("..//Data/rosaceae.csv") ## read in flower size, phenology and fruiot size data
hyster<-filter(hyster,genus=="Rhododendron")

colnames(hyster)[3]<-"specificEpithet"
plot(mean.rhody)
points(lonpoints,latpoints, col=c("royal blue","black"),pch=c(8,5),cex=0.7)
legend(-80,20,c("hysteranthous","seranthous"),pch=c(8,5),col=c("royal blue","black"))

rhody.data.2<-dplyr::left_join(rhody.data,hyster,by="specificEpithet")


d<-dplyr::filter(rhody.data.2,!is.na(hysteranthy)) ### 7668 rows have county coordiates


ggplot(d,aes(as.factor(hysteranthy),pdsi))+geom_boxplot()+theme_minimal()

