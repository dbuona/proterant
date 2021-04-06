#### Pull fruiting data for prunus
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library(tibble)
library(lubridate)
library(stringr)
library("ncdf4")
library(raster)
library(purrr)
library(tidyr)


setwd("~/Documents/git/proterant/investment/input")

d<-read.csv("midwest_round1Dec11.csv") # active datasheet

##subset to useful additions
d.add<-read.csv("species_additions - species_additions.csv")
d.add<-dplyr::filter(d.add,fruit=="Y")
unique(d.add$bbch.f)


##subset to flowering data
d.fruit<-dplyr::filter(d,fruit=="Y")
d.fruit<-rbind(d.fruit,d.add)
table(d.fruit$specificEpithet)

nested_prun <- d.fruit %>%
  group_by(specificEpithet)%>%   # prep for work by Species
  nest() %>%              # --> one row per Species
  ungroup() %>% 
  mutate(n = c(30,30,30,15,30,30,30,2,30,29,20,16,
               24))
sampled_pruny <- nested_prun%>%
  mutate(samp = map2(data, n, sample_n))

sampled_pruny<-sampled_pruny %>% 
  dplyr::select(-data) %>%
  unnest(samp)

sampled_pruny2<-dplyr::select(sampled_pruny,specificEpithet,id,references)

write.csv(sampled_pruny2,"fruit_measures.csv", row.names = FALSE)
