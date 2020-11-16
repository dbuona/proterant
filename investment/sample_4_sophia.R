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
d.geo<-filter(d,!is.na(lat))


nested_prun <- d.geo %>%
  group_by(specificEpithet) %>%   # prep for work by Species
  nest() %>%              # --> one row per Species
  ungroup() %>% 
  mutate(n = c(200,200,200,30,200,200,200,39,200,200,200,114,
               200,121,1))

print(nested_prun)
sampled_pruny <- nested_prun%>%
  mutate(samp = map2(data, n, sample_n))

sampled_pruny<-sampled_pruny %>% 
  dplyr::select(-data) %>%
  unnest(samp)

##remove duplicates
doney<-read.csv("pruno_checker.csv")

sampled_pruny<-filter(sampled_pruny,!references %in% doney$references)


write.csv(sampled_pruny,"midwest_round1.csv", row.names = FALSE)





