###A script to see if hysteranthy is changing in pep data


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


d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
unique(d$species)
### offset functional
d$offset<-d$l75.jd-d$fopn.jd
###0ffset physiological
d$offset2<-d$bb.jd-d$fbb.jd
###intermediate
d$offset3<-d$bb.jd-d$fopn.jd

sum<-d%>% group_by(species) %>% summarise(avg.offset = mean(offset,na.rm=TRUE))
sum2<-d%>% group_by(species) %>% summarise(avg.offset2 = mean(offset2,na.rm=TRUE))
sum3<-d%>% group_by(species) %>% summarise(avg.offset3 = mean(offset3,na.rm=TRUE))

a<-filter(d, species %in% c("ACRU","ACSA","AMSP" ,"BEAL" ,"BELE" ,"BEPA", "BEPO","FRAM","POTR" ,"QURU" ,"QUVE"))

#unique(hys$species) ##15 hysteranthous species (pro and syn)



ggplot(a,aes(year,offset))+geom_point()+facet_wrap(~species)
a1<-filter(a, species %in% c("ACPE","ACRU","ACSA","BEAL" ,"FRAM" ,"QURU"))
ggplot(a1,aes(year,offset))+geom_smooth(method="lm",aes(group=tree.id,col=species),linetype="dashed",se=FALSE,size=0.5)+geom_smooth(method="lm",aes(color=species),size=1.5,se=FALSE)+geom_smooth(method="lm",color="red",size=2,linetype="twodash")+theme_tufte()+scale_color_brewer(palette="Dark2")
ggplot(a1,aes(year,offset2))+geom_smooth(method="lm",aes(group=tree.id,col=species),linetype="dashed",se=FALSE,size=0.5)+geom_smooth(method="lm",aes(color=species),size=1.5,se=FALSE)+geom_smooth(method="lm",color="red",size=2,linetype="twodash")+theme_tufte()+scale_color_brewer(palette="Dark2")




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
fullfrax<-compy
fullfrax$taxa<-"Fraxinus excelsior"

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
fullaes<-compy
fullaes$taxa<-"Aesculus hippocastanum"
compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)

aesagg <- aggregate(compy[("year")], compy[c("s_id", "lat", "lon", "alt")],
                     FUN=length)


aes50 <- subset(aesagg, year>50)
compy<- compy[which(compy$s_id %in% aes50$s_id),]



write.csv(compy,"aes_delta_hyst.csv")

plotty<-ggplot(compy,aes(as.numeric(year),offset))+geom_smooth(method="lm",aes(group=locale),se=FALSE,color="black", linetype="dotted")+geom_smooth(method="lm",color="red")+ggtitle("Aesculus")


####Alnus glutinosa
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
fullalnus<-compy
fullalnus$taxa<-"Alnus glutinosa"

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)
alnusagg <- aggregate(compy[("year")], compy[c("s_id", "lat", "lon", "alt")],
                     FUN=length)


alnus50 <- subset(alnusagg, year>50)
compy<- compy[which(compy$s_id %in% alnus50$s_id),]


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
fullfag<-compy
fullfag$taxa<-"Fagus sylvatica"

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)
fagagg <- aggregate(compy[("year")], compy[c("s_id", "lat", "lon", "alt")],
                      FUN=length)


fag50 <- subset(fagagg, year>50)
compy<- compy[which(compy$s_id %in% fag50$s_id),]


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

fulltil<-compy
fulltil$taxa<-"Tilia cordata"

compy<-filter(compy,offset<70)
compy<-filter(compy,offset>(-70))
table(compy$year)
tiliaagg <- aggregate(compy[("year")], compy[c("s_id", "lat", "lon", "alt")],
                    FUN=length)


til50 <- subset(tiliaagg, year>50)
compy<- compy[which(compy$s_id %in% til50$s_id),]

write.csv(compy,"tilia_delta_hyst.csv")

###Corylus avenula
cor1<-read.csv("buo_107_000_011.csv", header=TRUE)
cor2<-read.csv("buo_107_000_060.csv",header=TRUE)
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
fullcory<-compy
fullcory$taxa<-"Corylus avenula"


##Betula pendula
cor1<-read.csv("buo_106_020_011.csv", header=TRUE)
cor2<-read.csv("buo_106_020_060.csv",header=TRUE)
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
fullbet<-compy
fullbet$taxa<-"Betula pendula"


###Acer platanus
cor1<-read.csv("buo_115_030_011.csv", header=TRUE)
cor2<-read.csv("buo_115_030_060.csv",header=TRUE)
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
fullacer<-compy
fullacer$taxa<-"Acer platinoides"

###Sorbus
cor1<-read.csv("buo_126_000_011.csv", header=TRUE)
cor2<-read.csv("buo_126_000_060.csv",header=TRUE)
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
fullsorb<-compy
fullsorb$taxa<-"Sorbus acuparia"

#Betula sp
cor1<-read.csv("buo_106_XXX_011.csv", header=TRUE)
cor2<-read.csv("buo_106_XXX_060.csv",header=TRUE)
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
fullbetquest<-compy
fullbetquest$taxa<-"Betula spp."

##syringa vulgaris
cor1<-read.csv("buo_127_000_011.csv", header=TRUE)
cor2<-read.csv("buo_129_000_060.csv",header=TRUE)
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
fullsyr<-compy
fullsyr$taxa<-"Syringa vulgaris"


###Tilia plataphyllos
cor1<-read.csv("buo_129_071_011.csv", header=TRUE)
cor2<-read.csv("buo_129_071_060.csv",header=TRUE)
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
fulltil.plat<-compy
fulltil.plat$taxa<-"Tilia Platyphyllos"


##Quercus spp.
cor1<-read.csv("buo_111_xxx_011.csv", header=TRUE)
cor2<-read.csv("buo_111_xxx_060.csv",header=TRUE)
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
fullquer<-compy
fullquer$taxa<-"Quercus spp."

###Prunis avium
cor1<-read.csv("buo_222_xxx_011.csv", header=TRUE)
cor2<-read.csv("buo_222_xxx_060.csv",header=TRUE)
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
fullpru<-compy
fullpru$taxa<-"Prunus avium"
big.dat<-rbind(fullacer,fullaes,fullalnus,fullbet,fullbetquest,fullcory,fullfag,fullfrax,fullpru,fullquer,fullsorb,fullsyr,fulltil,fulltil.plat)

###remove NA's
big.dat<-filter(big.dat, !is.na(offset))

##remove offsets greater that 90
big.dat<-filter(big.dat, offset<90)
big.dat<-filter(big.dat, offset>-90)

write.csv(big.dat,"PEP_OFFSET_FUll.csv", row.names = FALSE)


