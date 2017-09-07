
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Desktop/Usable Pep")

library(dplyr)
library(tidyr)
library(ggplot2)



#library(car)
#scatterplot(day~ lat| bbch, data=oneoone, main="Pep 101") 



cor1<-read.csv("buo_102_040_011.csv", header=TRUE)
cor2<-read.csv("buo_102_040_060.csv",header=TRUE)


names(cor1)[names(cor1)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor1<-separate(cor1,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(sep_rec1)

names(cor2)[names(cor2)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor2<-separate(cor2,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(cor2)

cor<-rbind(cor2,cor1)
table(cor$year)

cor<-filter(cor, year>2010)

ggplot(cor,aes(as.numeric(lat),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Alnus glutinosa across latitude")
