rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")
library(ggplot2)
library(tidyr)
library(dplyr)
library(chillR)
library(brms)
library(tibble)
library(lme4)
library("lmerTest")

HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
#HF2<-read.csv("HarvardForest/hf003-06-mean-spp.csv",header=TRUE)



##Hypothesis 1: variation in FLS is a product of variation in climate between phases
#https://cran.r-project.org/web/packages/chillR/vignettes/hourly_temperatures.html
weather<-read.csv("..//FLOBUDS/data/hf000-01-daily-m.csv",header = TRUE)
weather<-dplyr::select(weather,c("date","airtmax","airtmin"))
weather<-separate(weather,date,c("Year","Month","Day"),sep="-",remove=TRUE)
colnames(weather)<-c("Year","Month","Day","Tmax","Tmin")
sapply(weather,mode) #mode(weather)
weather$Year<-as.numeric(weather$Year)
weather$Month<-as.numeric(weather$Month)
weather$Day<-as.numeric(weather$Day)

unique(weather$Year)
weather<-filter(weather,Year>=1989)

chambercorr<-weather
chambercorr$Tmax<-4
chambercorr$Tmin<-4

unique(chambercorr$Year)
chambercorr<-filter(chambercorr,Year>=2002)

#all_daylengths<-cbind(JDay=1:365,sapply(daylength(latitude=42.5,JDay=1:365),cbind)) ## calculate day length at HF on every day of year
#ad<-as.data.frame(all_daylengths) ## data frame of day length


chambercorr<-make_all_day_table(chambercorr)

hourtemps<-stack_hourly_temps(chambercorr, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)



chillref<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=60))
chillrefshort<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=30))
###3 now compare to real data
unique(weather$Year)
weather<-filter(weather,Year>=1989)

#all_daylengths<-cbind(JDay=1:365,sapply(daylength(latitude=42.5,JDay=1:365),cbind)) ## calculate day length at HF on every day of year
#ad<-as.data.frame(all_daylengths) ## data frame of day length


weather<-make_all_day_table(weather)

hourtemps<-stack_hourly_temps(weather, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)

climatt<-as.data.frame(chilling(hourtemps,Start_JDay=365-31-30-15,End_JDay=31+28+31+15)) ##Oct 15-April 15

envtoreal<-data.frame(UtahEst=c(paste(mean(climatt$Utah_Model),((sd(climatt$Utah_Model))),sep=","),chillref$Utah_Model,chillrefshort$Utah_Model),
                      CpEst=c(paste(mean(climatt$Chill_portions),sd(cllimate),chillref$Chill_portions,chillrefshort$Chill_portions),
                      ChEst=c(mean(climatt$Chilling_Hours),chillref$Chilling_Hours,chillrefshort$Chilling_Hours),
           treatment=c("HFaverage Oct 15-April 15","60 days chamber","30 days chamber")))


able<-data.frame(Estimate=c("Utah Model","Chill Hours", "Dynamic Model"),
Havard_Forest=c(mean(climatt$Utah_Model),mean(climatt$Chilling_Hours),mean(climatt$Chill_portions)),
sd=c(sd(climatt$Utah_Model),sd(climatt$Chilling_Hours),sd(climatt$Chill_portions)),
Chamber_30_days=c(chillrefshort$Utah_Model,chillrefshort$Chilling_Hours,chillrefshort$Chill_portions),
Chamber_60_days=c(chillref$Utah_Model,chillref$Chilling_Hours,chillref$Chill_portions))

is.num <- sapply(able, is.numeric)
able[is.num] <- lapply(able[is.num], round, 2)
save.image("expfieldcomp.Rda")
write.csv(able,"..//FLOBUDS/input/expfieldcomparisontable.csv",row.names = TRUE)

           

