###mergine two days of field sampling for flobuds project
### Dan B on Tuesday Aug 8

rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(dplyr)

setwd("~/Documents/git/proterant/FLOBUDS")

d1<-read.csv("DanB_HFPOINTS_JUL_2017.csv", header=TRUE)
d2<-read.csv("DB_HF_Points_AUG_2017.csv",header=TRUE)

fulldat<-full_join(d2,d1, by="name")
fulldat<-select(fulldat, name,lat.x,lon.x,lat.y,lon.y,ele.x,DBH.y,locality.y,CreationTime.x,CreatiionTime.y)

write.csv(fulldat,"flobud_datapoints.csv",row.names=FALSE)
### after this step, new dbh's were added directly in excel, so you should not need to reproduce this code ever