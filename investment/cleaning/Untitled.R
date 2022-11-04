###make FNA.prunus

### Started by Dan 07 Oct 2002
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
options(mc.cores = parallel::detectCores())


library(dplyr)
library("housingData")
library(stringr)
library("ncdf4")
library(raster)
library(ggrepel)
library(patchwork)
library(brms)
library("grid")
library(tidybayes)
library(bayesplot)
library("caper")
library(phytools)
library(geiger)


#-------------------------#
#------Part 1--------------#
#--------------------------#
setwd("~/Documents/git/proterant/investment")
FNA<-read.csv("Data/cherry_data.csv") ## measurement and FLS data from FNA
##clean specicies name

FNA$species<-ifelse(FNA$species=="hortunlana","hortulana",FNA$species)
FNA$species<-ifelse(FNA$species=="gladulosa","glandulosa",FNA$species)
FNA$species<-ifelse(FNA$species=="speciosa","speciosa",FNA$species)
FNA$species<-ifelse(FNA$species=="fasiculata","fasciculata",FNA$species)

#### get herbaria specimen coordinates
mid.herb<-read.csv("SymbOutput_2020-10-26_133412_DwC-A/occurrences.csv")
mid.herb<-filter(mid.herb,specificEpithet %in% unique(FNA$species)) ## filter to Species we have data 4
unique(mid.herb$specificEpithet)
## select rows that already have coordinate
pruno.ref<-filter(mid.herb,!is.na(decimalLatitude))
table(pruno.ref$stateProvince)
table(pruno.ref$specificEpithet)
pruno.ref$lon<-pruno.ref$decimalLongitude
pruno.ref$lat<-pruno.ref$decimalLatitude

##now add county level coordiantes
#### two without coordinates
pruno.unref<-filter(mid.herb,is.na(decimalLatitude)) ## filter entries with no lat/lon

### prep county data
geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" 

colnames(geoCounty)[c(2,7)]<-c("County_name","county") ## MATCH NAMES

pruno.unref<-filter(pruno.unref,id!=17453277) #remvoe problematic canadian entry
pruno.unref$county<-tolower(pruno.unref$county) ## make ours lower case



pruno.unref$county<-gsub("county","",pruno.unref$county) ### get rid of extenious coounty info
pruno.unref$county<-gsub("()","",pruno.unref$county) # ""
pruno.unref$county<-gsub("co.","",pruno.unref$county,fixed = TRUE) # ""




pruno.unref<-dplyr::left_join(pruno.unref,geoCounty,by=c("county","stateProvince"))
pruno.unref<-dplyr::select(pruno.unref,-County_name,-fips, -state)
intersect(colnames(pruno.unref),colnames(pruno.ref))

pruno.ref$geoclass<-"geo_referenced"
pruno.unref$geoclass<-"county_centroid"
pruno.unref<-filter(pruno.unref,!is.na(lon)) # select only entries with lat lon
pruneo<-rbind(pruno.ref,pruno.unref)

####calculated mean pdsi for each species
palmer.b <- brick("Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-pruneo$lon # make vector of prunus coordinates
latpoints<-pruneo$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
sd.prunus <-calc(palmer.b, fun = sd,na.rm=TRUE) #
min.prunus <-calc(palmer.b, fun = min,na.rm=TRUE) 
ext<-raster::extract(mean.prunus,extract.pts,method="simple")
ext2<-raster::extract(sd.prunus,extract.pts,method="simple")
ext3<-raster::extract(min.prunus,extract.pts,method="simple")
pruneo$pdsi<-ext
pruneo$pdsi.sd<-ext2
pruneo$pdsi.min<-ext3


pdsi<- pruneo %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanpdsi=mean(pdsi,na.rm=TRUE),sdmean=sd(pdsi,na.rm=TRUE),minpdsi=mean(pdsi.min,na.rm=TRUE),sdmin=sd(pdsi.min,na.rm=TRUE))
colnames(pdsi)[1]<-"species"
setdiff(FNA$species,pdsi$species) ## lose speciosa


## combine data
FNA<-left_join(FNA,pdsi)
##quick plots

FNA$FLSnum<-NA
FNA$FLSnum[FNA$FLS=="before"]<-1
FNA$FLSnum[FNA$FLS=="before/with"]<-2
FNA$FLSnum[FNA$FLS=="with"]<-3
FNA$FLSnum[FNA$FLS=="after"]<-4
FNA$FLSnum<-as.integer(FNA$FLSnum)

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
FNA$meanpdsi.z<-zscore(FNA$meanpdsi)
FNA$minpdsi.z<-zscore(FNA$minpdsi)
FNA$petal.z<-zscore(FNA$petal_high)
FNA$inflor.z<-zscore(FNA$inflor_high)
FNA$fruit.z<-zscore(FNA$fruit_high)

FNA$logFLS<-log(FNA$FLSnum)

write.csv(FNA,"Input/input_clean/FNA_final.csv")
