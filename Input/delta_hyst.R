###A script to see if hysteranthy is changing in pep data

rm(list=ls()) 
options(stringsAsFactors = FALSE)
getOption("device")
options(device="RStudioGD")
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

d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
unique(d$species)

a<-filter(d, species %in% c("ACRU","ACSA","AMSP" ,"BEAL" ,"BELE" ,"BEPA", "BEPO","FRAM","POTR" ,"QURU" ,"QUVE"))

#unique(hys$species) ##15 hysteranthous species (pro and syn)
a$offset<-a$l75.jd-a$fopn.jd
###budburst to 
a$offset2<-a$bb.jd-a$fopn.jd


ggplot(a,aes(year,offset))+geom_point()+facet_wrap(~species)
a1<-filter(a, species %in% c("ACPE","ACRU","ACSA","BEAL" ,"FRAM" ,"QURU"))
ggplot(a1,aes(year,offset))+geom_smooth(method="lm",aes(group=species,linetype=species),color="black")+geom_smooth(method="lm",color="red",size=1)+theme_base()
ggplot(a1,aes(year,offset2))+geom_smooth(method="lm",aes(group=tree.id,col=species),se=FALSE,size=0.25)+geom_smooth(method="lm",aes(color=species),size=1.5,se=FALSE)+geom_smooth(method="lm",color="red",size=2,linetype="dashed")+theme_base()+scale_color_brewer(palette="Dark2")

#a<-gather(a, phase,day,4:7)
#a<-filter(a, phase %in% c("fopn.jd","l75.jd"))
#bet<-filter(a, species %in% c("BEAL","BEPA","BEPO"))
#ggplot(bet,aes(year,day))+geom_point(aes(shape=phase))+geom_smooth(method='lm',aes(,color=phase))+theme_base()+facet_wrap(~species)
#bet$lat<-"42.5315"
#bet$lon<-"72.19"
#bet<-unite(bet,locale,lat,lon,sep=",",remove=FALSE)
bet$bbch<-ifelse(bet$phase=="fopn.jd","60","11")
bet$year<-as.character(bet$year)
bet$day<-as.character(bet$day)
bet$s_id<-bet$species
#############################Alnus glutinosa trends in pep

###Fraxinus
cor1<-read.csv("buo_120_000_011.csv", header=TRUE)
cor2<-read.csv("buo_120_000_060.csv",header=TRUE)
#cor3<-read.csv("buo_106_020_010.csv",header=TRUE)

names(cor1)[names(cor1)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor1<-separate(cor1,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(sep_rec1)

names(cor2)[names(cor2)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor2<-separate(cor2,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(cor2)


cor<-rbind(cor2,cor1)

cor<-unite(cor,locale,lat,lon,sep=",",remove=FALSE)

compy<-cor
#ggplot(compy,aes(as.numeric(year),as.numeric(day)))+geom_point(aes())+geom_smooth(method='lm',aes(col=bbch))+theme_base()+ggtitle("Betula phenology since since 1960")+facet_wrap(~locale)
compy$pheno<-ifelse(compy$bbch=="11","leaf","flower")
compy<-dplyr::select(compy,-bbch)
compy<-spread(compy,pheno,day)
compy$leaf<-as.numeric(compy$leaf)
compy$flower<-as.numeric(compy$flower)

compy$offset<-compy$leaf-compy$flower

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)
fraxagg <- aggregate(compy[("year")], compy[c("s_id", "lat", "lon", "alt")],
                    FUN=length)


frax50 <- subset(fraxagg, year>50)
compy<- compy[which(compy$s_id %in% frax50$s_id),]




#compy<-filter(compy, s_id %in% c(range)) ##stations have data 1960-2000
#compy<-filter(compy, year>1960)
write.csv(compy,"fraxinus_delta_hyst.csv")

plotty<-ggplot(compy,aes(as.numeric(year),offset))+geom_smooth(method="lm",aes(group=locale),se=FALSE,color="black", linetype="dotted")+geom_smooth(method="lm",color="red")+ggtitle("Fraxinus")

###########Aesculus

cor1<-read.csv("buo_101_000_011.csv", header=TRUE)
cor2<-read.csv("buo_101_000_060.csv",header=TRUE)
#cor3<-read.csv("buo_106_020_010.csv",header=TRUE)

names(cor1)[names(cor1)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor1<-separate(cor1,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(sep_rec1)

names(cor2)[names(cor2)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor2<-separate(cor2,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(cor2)


cor<-rbind(cor2,cor1)

cor<-unite(cor,locale,lat,lon,sep=",",remove=FALSE)

compy<-cor
#ggplot(compy,aes(as.numeric(year),as.numeric(day)))+geom_point(aes())+geom_smooth(method='lm',aes(col=bbch))+theme_base()+ggtitle("Betula phenology since since 1960")+facet_wrap(~locale)
compy$pheno<-ifelse(compy$bbch=="11","leaf","flower")
compy<-dplyr::select(compy,-bbch)
compy<-spread(compy,pheno,day)
compy$leaf<-as.numeric(compy$leaf)
compy$flower<-as.numeric(compy$flower)

compy$offset<-compy$leaf-compy$flower

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)

aesagg <- aggregate(compy[("year")], compy[c("s_id", "lat", "lon", "alt")],
                     FUN=length)


aes50 <- subset(aesagg, year>50)
compy<- compy[which(compy$s_id %in% aes50$s_id),]



write.csv(compy,"aes_delta_hyst.csv")

plotty<-ggplot(compy,aes(as.numeric(year),offset))+geom_smooth(method="lm",aes(group=locale),se=FALSE,color="black", linetype="dotted")+geom_smooth(method="lm",color="red")+ggtitle("Aesculus")

cor1<-read.csv("buo_102_040_011.csv", header=TRUE)
cor2<-read.csv("buo_102_040_060.csv",header=TRUE)
#cor3<-read.csv("buo_106_020_010.csv",header=TRUE)

names(cor1)[names(cor1)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor1<-separate(cor1,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(sep_rec1)

names(cor2)[names(cor2)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor2<-separate(cor2,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(cor2)


cor<-rbind(cor2,cor1)

cor<-unite(cor,locale,lat,lon,sep=",",remove=FALSE)

compy<-cor
#ggplot(compy,aes(as.numeric(year),as.numeric(day)))+geom_point(aes())+geom_smooth(method='lm',aes(col=bbch))+theme_base()+ggtitle("Betula phenology since since 1960")+facet_wrap(~locale)
compy$pheno<-ifelse(compy$bbch=="11","leaf","flower")
compy<-dplyr::select(compy,-bbch)
compy<-spread(compy,pheno,day)
compy$leaf<-as.numeric(compy$leaf)
compy$flower<-as.numeric(compy$flower)

compy$offset<-compy$leaf-compy$flower

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))

###select only stations with data from 1961 and 2000 which i dont need I think if fitting mixed model
sixty<-filter(compy, year==1951)
oh<-filter(compy,year==2000)


sixty<-filter(sixty, !is.na(offset))
sitystations<-sixty$s_id
Ohstations<-oh$s_id
range<-intersect(sitystations,Ohstations)

compy<-filter(compy, s_id %in% c(range)) ##stations have data 1960-2000
compy<-filter(compy, year>1960)

write.csv(compy,"alnus_delta_hyst.csv")

ggplot(compy,aes(as.numeric(year),offset))+geom_smooth(method="lm",aes(group=s_id),se=FALSE,color="black", linetype="dotted")+geom_smooth(method="lm",color="red")+ggtitle("Alnus glutinosa")+theme_base()


#####fagus
cor1<-read.csv("buo_108_010_011.csv", header=TRUE)
cor2<-read.csv("buo_108_010_060.csv",header=TRUE)
#cor3<-read.csv("buo_106_020_010.csv",header=TRUE)

names(cor1)[names(cor1)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor1<-separate(cor1,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(sep_rec1)

names(cor2)[names(cor2)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor2<-separate(cor2,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(cor2)


cor<-rbind(cor2,cor1)

cor<-unite(cor,locale,lat,lon,sep=",",remove=FALSE)

compy<-cor
#ggplot(compy,aes(as.numeric(year),as.numeric(day)))+geom_point(aes())+geom_smooth(method='lm',aes(col=bbch))+theme_base()+ggtitle("Betula phenology since since 1960")+facet_wrap(~locale)
compy$pheno<-ifelse(compy$bbch=="11","leaf","flower")
compy<-dplyr::select(compy,-bbch)
compy<-spread(compy,pheno,day)
compy$leaf<-as.numeric(compy$leaf)
compy$flower<-as.numeric(compy$flower)

compy$offset<-compy$leaf-compy$flower

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)
seven<-filter(compy, year==1951)
oh<-filter(compy,year==1990) #### not enough for 2000


seven<-filter(seven, !is.na(offset))
sitystations<-seven$s_id
Ohstations<-oh$s_id
range<-intersect(sitystations,Ohstations)

compy<-filter(compy, s_id %in% c(range)) ##stations have data 1960-2000
compy<-filter(compy, year>1960)

write.csv(compy,"fag_delta_hyst.csv")

#####Tilia
cor1<-read.csv("buo_129_070_011.csv", header=TRUE)
cor2<-read.csv("buo_129_070_060.csv",header=TRUE)
#cor3<-read.csv("buo_106_020_010.csv",header=TRUE)

names(cor1)[names(cor1)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor1<-separate(cor1,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(sep_rec1)

names(cor2)[names(cor2)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #alternative method
#try seperation again
cor2<-separate(cor2,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" )
head(cor2)


cor<-rbind(cor2,cor1)

cor<-unite(cor,locale,lat,lon,sep=",",remove=FALSE)

compy<-cor
#ggplot(compy,aes(as.numeric(year),as.numeric(day)))+geom_point(aes())+geom_smooth(method='lm',aes(col=bbch))+theme_base()+ggtitle("Betula phenology since since 1960")+facet_wrap(~locale)
compy$pheno<-ifelse(compy$bbch=="11","leaf","flower")
compy<-dplyr::select(compy,-bbch)
compy<-spread(compy,pheno,day)
compy$leaf<-as.numeric(compy$leaf)
compy$flower<-as.numeric(compy$flower)

compy$offset<-compy$leaf-compy$flower

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)
seven<-filter(compy, year==1951)
oh<-filter(compy,year==1990) #### not enough for 2000


seven<-filter(seven, !is.na(offset))
sitystations<-seven$s_id
Ohstations<-oh$s_id
range<-intersect(sitystations,Ohstations)

compy<-filter(compy, s_id %in% c(range)) ##stations have data 1960-2000
compy<-filter(compy, year>1960)

write.csv(compy,"tilia_delta_hyst.csv")

