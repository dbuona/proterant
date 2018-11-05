###PEP analy sis, use pep to find trends in hysteranthy
rm(list=ls()) 
options(stringsAsFactors = FALSE)
getOption("device")
options(device="quartz") ###if at arb if not RStudioGD
sessionInfo()

setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)
library(ggthemes)
library(brms)
big.dat<-read.csv("PEP_OFFSET_FULL.csv")
#big.dat<-filter(big.dat,offset<42)
#big.dat<-filter(big.dat,offset>(-42))

###choose stations with at least 10 years of data
sagg<-aggregate(big.dat[("year")], big.dat[c("s_id", "lat", "lon", "alt")],
                    FUN=length)


big.dat10 <- subset(sagg, year>5)
big.dat<- big.dat[which(big.dat$s_id %in% big.dat10$s_id),]


sum<-big.dat%>% group_by(s_id,taxa,lat,alt) %>% summarise(avg.offset = mean(offset))
?summarise_at()
###does latitude or altitude predict offset
ggplot(sum,aes(lat,avg.offset))+geom_point(aes(color=taxa))+geom_hline(yintercept=0)+geom_smooth(method="lm")+facet_wrap(~taxa)
ggplot(sum,aes(alt,avg.offset))+geom_point(aes(color=taxa))+geom_hline(yintercept=0)+geom_smooth(method="lm")+facet_wrap(~taxa)



priorz<-get_prior(offset~alt+lat, data=big.dat)
cline.brms<-brm(offset~alt+lat+(alt+lat|taxa),data=big.dat, cores=4)
               

simple.cline<-brm(offset~alt+lat,data=big.dat)
