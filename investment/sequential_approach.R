#### simple models, no jointiness
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


setwd("~/Documents/git/proterant/investment/input")
load("sequentialmodeling.Rda")
d<-read.csv("midwest_round1Dec11.csv") # active datasheet

##subset to useful additions
d.add<-read.csv("species_additions - species_additions.csv")
d.add<-dplyr::filter(d.add,flowers=="Y")
unique(d.add$bbch.f)
##subset to flowering data
d.flo<-dplyr::filter(d,flowers=="Y")

###cleanup #clean up everybody
d.flo<-dplyr::filter(d.flo,bbch.f<80) #clean
d.flo<-dplyr::filter(d.flo,bbch.f>59) #clean
table(d.flo$specificEpithet)
table(d.add$specificEpithet)
d.flo<-rbind(d.flo,d.add)

d.flo$bbch.v<-ifelse(is.na(d.flo$bbch.v),0,d.flo$bbch.v)


###from all the stuff below we know a rescaled model with doy centered as a predictor is best
####fuctions
extract_coefsF<-function(x){rownames_to_column(as.data.frame(coef(x, summary=TRUE,probs=c(0.025,0.25,0.75,0.975))),"species")
}


### now control for date of observation as a covariate
#if month was missing we impute april
d.flo$month<-ifelse(d.flo$month==0,4,d.flo$month)
d.flo$month<-ifelse(is.na(d.flo$month),4,d.flo$month)

#if day was missing we imput 15
d.flo$day<-ifelse(d.flo$day==0,15,d.flo$day)
d.flo$day<-ifelse(is.na(d.flo$day),15,d.flo$day)

d.flo$year<-ifelse(d.flo$day>32,d.flo$day,d.flo$year)
d.flo$day<-ifelse(d.flo$day>32,15,d.flo$day)

#clean year and impute 1980 if missing
d.flo$year<-ifelse(d.flo$year==2,2002,d.flo$year)
d.flo$year<-ifelse(d.flo$year==5,2005,d.flo$year)
d.flo$year<-ifelse(d.flo$year==6,2006,d.flo$year)
d.flo$year<-ifelse(d.flo$year==198,1980,d.flo$year)
d.flo$year<-ifelse(is.na(d.flo$year),1980,d.flo$year)
d.flo$year<-ifelse(d.flo$year==0,1980,d.flo$year)
unique(d.flo$year)

d.flo$eventDate2<-paste(d.flo$year,d.flo$month,d.flo$day,sep="-")
d.flo<-filter(d.flo,!is.na(eventDate2))
d.flo$eventDate3<-as.Date(d.flo$eventDate2,format =  "%Y-%m-%d")
d.flo$doy<-yday(d.flo$eventDate3)

d.flo$doy.cent<-d.flo$doy-mean(d.flo$doy,na.rm=TRUE)
d.flo<-filter(d.flo,!is.na(doy.cent))
table(d.flo$specificEpithet)
## now rescale BBCH so its even 
unique(d.flo$bbch.v)
d.flo$bbch.v.scale<-NA
d.flo$bbch.v.scale[d.flo$bbch.v==0]<-0
d.flo$bbch.v.scale[d.flo$bbch.v==7]<-1
d.flo$bbch.v.scale[d.flo$bbch.v==9]<-1
d.flo$bbch.v.scale[d.flo$bbch.v==10]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==11]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==14]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==15]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==17]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==19]<-5

###add pdsi
palmer.b <- brick("..//lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d$lon # make vector of prunus coordinates
latpoints<-d$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
sd.prunus <-calc(palmer.b, fun = sd,na.rm=TRUE) #
min.prunus <-calc(palmer.b, fun = min,na.rm=TRUE) 
ext<-raster::extract(mean.prunus,extract.pts,method="simple")
ext2<-raster::extract(sd.prunus,extract.pts,method="simple")
ext3<-raster::extract(min.prunus,extract.pts,method="simple")
d$pdsi<-ext
d$pdsi.sd<-ext2
d$pdsi.min<-ext3

### model control for doy,
mod_doy_gaus<-brm(bbch.v.scale~doy.cent+(1|specificEpithet),data=d.flo)

gaus.out<-extract_coefsF(mod_doy_gaus)
gaus.out<-dplyr::select(gaus.out,1:7)

colnames(gaus.out)<-c("specificEpithet","EstimatedFLS","FLSSE","Q2.5","Q25","Q75","Q97.5")
gaus.out$specificEpithet <- factor(gaus.out$specificEpithet, levels = gaus.out$specificEpithet[order(gaus.out$EstimatedFLS)])

d<-left_join(d,gaus.out)
mod_pdsi<-brm( EstimatedFLS | mi(FLSSE)~ pdsi,data=d)
summary(mod_pdsi)
fixef(mod_pdsi,probs = c(.025,.25,.75,.975))

mod_pdsi_min<-brm(pdsi.min~me(EstimatedFLS, FLSSE)+(1|specificEpithet),data=d)
conditional_effects(mod_pdsi,probs = c(.25,.75))
fixef(mod_pdsi,probs = c(.025,.25,.75,.975))

#mod_pdsi_min_sd<-brm(pdsi|mi(pdsi.sd)~me(EstimatedFLS, FLSSE)+(1|specificEpithet),data=d.flo)


flowers<-read.csv("flowermeasures - measurements.csv")
flower<-left_join(flowers,gaus.out)
class(flower$id)
mod_size_min<-brm(pental_lengh_mm~me(EstimatedFLS, FLSSE)+(1|id)+(1|specificEpithet),data=flower, iter=8000,warmup=7000, control=list(adapt_delta=.99))
fixef(mod_size_min,probs = c(.025,.25,.75,.975))
coef(mod_size_min,probs = c(.025,.25,.75,.975))
summary(mod_size_min)
save.image("sequentialmodeling.Rda")
