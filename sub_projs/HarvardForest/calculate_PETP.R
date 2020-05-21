rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_projs/")

library(ape)
library(phytools)
library(brms)
library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")
library("jpeg")
library("phylolm")
library(ggstance)
options(scipen = 999)
fluxy<-read.csv("HarvardForest/AMF_US-Ha2_BASE_HH_3-5.csv")
fluxy2<-read.csv("HarvardForest/AMF_US-Ha1_BASE_HR_14-5.csv")

colnames(fluxy)

#colnames(fluxy2)
fluxy2<-dplyr::select(fluxy2,TIMESTAMP_START,TIMESTAMP_END,P)
#head(fluxy2)
##flux data starts 2004-2019

fluxy<-dplyr::select(fluxy,TIMESTAMP_START,TIMESTAMP_END,LE_1_1_1,TA_1_1_1)

fluxy<-filter(fluxy,LE_1_1_1!=-9999)
head(fluxy$TIMESTAMP_START)
head(fluxy$TIMESTAMP_END)

library("bigleaf")
fluxy$ETP1<-LE.to.ET(c(fluxy$LE_1_1_1),c(fluxy$TA_1_1_1))
fluxy$year<- (fluxy$TIMESTAMP_START %/% 1e8)


fluxer<- fluxy %>% group_by(year) %>% summarise(meanETP=mean(ETP1))


#ETP1 is in kg m-2 s-1
##1 kg/m2/s = 86400 mm/day 
fluxer$ETP2<-fluxer$meanETP*86400*365 ### units mm/day#fluxy2<-read.csv("HarvardForest/AMF_US-Ha1_BASE_HR_14-5.csv")

###now get precip

#colnames(fluxy2)
fluxy2<-dplyr::select(fluxy2,TIMESTAMP_START,TIMESTAMP_END,P)
fluxy2<-filter(fluxy2,P!=-9999)

x <- 1293828893
fluxy2$year<- (fluxy2$TIMESTAMP_START %/% 1e8)

fluxer2<- fluxy2 %>% group_by(year) %>% summarise(TotalP=sum(P))
joint<-left_join(fluxer,fluxer2)

joint$PTEP<-joint$TotalP-joint$ETP2
joint$aridindex<-joint$ETP2/joint$TotalP

write.csv(joint,"PETP.HF.csv",row.names = FALSE)

