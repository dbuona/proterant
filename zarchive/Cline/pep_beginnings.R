
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Cline")
library(dplyr)
library(tidyr)
library(ggplot2)



#library(car)
#scatterplot(day~ lat| bbch, data=oneoone, main="Pep 101") 


##Alnus glutinosa
### alnus
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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>100)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Alnus glutinosa phenology since 1900")

####other species
cor1<-read.csv("buo_107_000_011.csv", header=TRUE)
cor2<-read.csv("buo_107_000_060.csv",header=TRUE)


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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>50)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Corylus avenula phenology since 1900")

####other species

####other species
cor1<-read.csv("buo_120_000_011.csv", header=TRUE)
cor2<-read.csv("buo_120_000_060.csv",header=TRUE)


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


colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>128)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(4495))
filted2<-filter(cor, s_id %in% c(4172))
table(filted$bbch)

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Fraxinus excelsior phenology since 1900")
ggplot(filted2,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Fraxinus excelsior phenology since 1900")

##dig more into faxinus
table(cor$s_id)


####other species
cor1<-read.csv("buo_106_xxx_011.csv", header=TRUE)
cor2<-read.csv("buo_106_xxx_060.csv",header=TRUE)


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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)
most<-filter(most,n>100)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Betula pendula phenology since 1961, Zagreb, Croatia")

####other species

cor1<-read.csv("buo_108_010_011.csv", header=TRUE)
cor2<-read.csv("buo_108_010_060.csv",header=TRUE)


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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>50)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Fagus phenology since 1900")

####other species

cor1<-read.csv("buo_107_000_011.csv", header=TRUE)
cor2<-read.csv("buo_107_000_060.csv",header=TRUE)


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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>50)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Corylus avenula phenology since 1900")

####other species

cor1<-read.csv("buo_107_000_011.csv", header=TRUE)
cor2<-read.csv("buo_107_000_060.csv",header=TRUE)


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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>50)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Corylus avenula phenology since 1900")

####other species

cor1<-read.csv("buo_107_000_011.csv", header=TRUE)
cor2<-read.csv("buo_107_000_060.csv",header=TRUE)


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

cor<-filter(cor, year>1900)
colnames(cor)
cor<-unite(cor,locale,lat,lon,sep=",")
most<-count(cor,s_id)

most<-filter(most,n>50)
stations<-most$s_id
filted<-filter(cor, s_id %in% c(stations))

#ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(linetype=bbch,color=s_id))
ggplot(filted,aes(as.numeric(year),as.numeric(day)))+geom_point(aes(shape=bbch))+geom_smooth(method='lm',aes(col=bbch))+theme()+ggtitle("Corylus avenula phenology since 1900")

####other species









