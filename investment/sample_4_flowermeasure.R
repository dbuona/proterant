### Explore the prunus data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(stringr)
library("ncdf4")
library(raster)
library(ggplot2)
library("brms")
library(dplyr)
library(purrr)
library(tidyr)
setwd("~/Documents/git/proterant/investment/input")

d<-read.csv("midwest_round1Dec11.csv")
d.measure<-dplyr::filter(d,bbch.f==65)
table(d.measure$specificEpithet)

nested_prun <- d.measure %>%
  group_by(specificEpithet) %>%   # prep for work by Species
  nest() %>%              # --> one row per Species
  ungroup() %>% 
  mutate(n = c(30,30,30,8,30,30,30,6,30,30,30,7,
               30))

print(nested_prun)
sampled_pruny <- nested_prun%>%
  mutate(samp = map2(data, n, sample_n))

sampled_pruny<-sampled_pruny %>% 
  dplyr::select(-data) %>%
  unnest(samp)

sampled_pruny2<-dplyr::select(sampled_pruny,specificEpithet,id,references)

write.csv(sampled_pruny2,"flowermeasures.csv", row.names = FALSE)
