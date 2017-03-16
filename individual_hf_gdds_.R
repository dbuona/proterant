###for indv species
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

setwd("~/Documents/git/proterant")
hf<- read.csv("data/WeatherData.csv", header=TRUE)
hf<-distinct(hf,Date,.keep_all = TRUE)
d<-read.csv("data/hf003-05-mean-ind.csv",header=TRUE)
###adding gdd
hf$gdd <- hf$AirT - 5
hf$gdd <-ifelse(hf$gdd>0, hf$gdd, 0)
hf$gdd <-ifelse(!is.na(hf$gdd), hf$gdd, 0) #added by dan-this is the problem line
hf$count <- ave(
  hf$gdd, hf$Year, 
  FUN=function(x) cumsum(c(0, head(x, -1)))
)
hf<-filter(hf,Year>=1990)
###mergine datasets
hf<-rename(hf, year=Year)
#d<-d%>%
  #dplyr::select(year,species, bb.jd,l75.jd,fbb.jd,fopn.jd)
d<- as.data.frame(rapply(object = d, f = round, classes = "numeric", how = "replace", digits = 0)) 
df<-left_join(hf, d)
##subsetting
df2<-df%>%
  dplyr::select(year,species, JD, bb.jd, count) %>%

df2$day<- ifelse(df2$JD==df2$bb.jd,df2$JD,NA)
df2<-na.omit(df2)
