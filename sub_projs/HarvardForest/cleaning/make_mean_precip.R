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

weather<-read.csv("..//FLOBUDS/data/hf000-01-daily-m.csv",header = TRUE)
weather<-dplyr::select(weather,c("date","prec"))
weather<-separate(weather,date,c("Year","Month","Day"),sep="-",remove=TRUE)

weather.means<- weather %>% group_by(Year) %>% summarise(AP=sum(prec,na.rm=TRUE))
write.csv(weather.means, "..//sub_projs/HarvardForest/mean.HF.precip.csv",row.names = FALSE)

