###map making for prunus manuscript updated by Dan on Jan 18
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library("rstan")

library(phytools)
library(ape)
library(lubridate)
library(stringr)
library("tidybayes")
library(raster)

require(mapdata); require(maptools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/Documents/git/proterant/investment/Input")

d<-read.csv("input_clean/pruno_clean_pdsi.csv")
d<-dplyr::filter(d,specificEpithet!="murrayana")
unique(d$specificEpithet)
###remove weird outlyers
d<-filter(d,lon<(-60))
d<-filter(d,lon>(-150))

usa <- map_data("usa")


jpeg("..//Plots/map.jpeg", width=11, height=7,unit="in",res=300)
ggplot()+geom_polygon(data=usa,aes(long,lat,group=group),fill="white",color="black")+
  geom_point(data=d,aes(lon,lat,color=specificEpithet,shape=specificEpithet),alpha=0.9)+ggthemes::theme_few()+
  scale_color_viridis_d(option="turbo")+
  scale_shape_manual(values = rep(c(0,1,2,15,16,17), len = 14))
dev.off()

