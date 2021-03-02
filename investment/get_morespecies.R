###this code add the ungeoreferences sheets for the 3 species I am lacking in good quantity.

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
library(dplyr)
library(purrr)
library(tidyr)

d<-read.csv("Data/midwestherbaria_data_sheet.csv")

set.seed(100)
##start with georeferences
d.nogeo<-filter(d,is.na(lat))
d.need<-filter(d.nogeo,specificEpithet %in% c("texana","subcordata","alleghaniensis"))


write.csv(d.need,"species_additions.csv", row.names = FALSE)
