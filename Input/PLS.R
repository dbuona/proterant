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
library(RColorBrewer)
library(ggstance)
set.seed(613)
HF<-read.csv("HarvardForest/hf003-06-mean-spp.csv",header=TRUE)
#HF2<-read.csv("HarvardForest/hf003-06-mean-spp.csv",header=TRUE)
#load("fieldexamples")



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
weather<-filter(weather,Year>=1990)

weather<-make_all_day_table(weather)


Acrb<-filter(HF, species=="ACRU")
Acrb<-filter(Acrb, year<=2002)

burst<-dplyr::select(Acrb,year,bb.jd)
colnames(burst)<-c("Year","pheno")
flo<-dplyr::select(Acrb,year,fopn.jd)
colnames(flo)<-c("Year","pheno")


weather<-fix_weather(weather[which(weather$Year>1990),])
#Plots look much better with weather<-fix_weather(KA_weather)
dc<-daily_chill(stack_hourly_temps(weather,50.4), 11)

plscf.bb<-PLS_chill_force(daily_chill_obj=dc, bio_data_frame=burst, split_month=5,end_at_pheno_end=TRUE,chill_models=c("Utah_Chill_Units"))

plscf.flo<-PLS_chill_force(daily_chill_obj=dc, bio_data_frame=flo, split_month=5,end_at_pheno_end=TRUE,chill_models=c("Utah_Chill_Units"))







PLS_results_path<-paste(getwd(),"/PLS_output.bb",sep="")
PLS_results_path2<-paste(getwd(),"/PLS_output.flo",sep="")

plot_PLS(plscf.bb,PLS_results_path)

plot_PLS(plscf.flo,PLS_results_path2)
