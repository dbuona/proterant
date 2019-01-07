rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(tidyverse)
library(rworldmap)
library(lme4)
library(car)
library(rgdal)


#read in clean data of offset for each species
aln<-read.csv("alnus_delta_hyst.csv",header=TRUE)
frax<-read.csv("fraxinus_delta_hyst.csv",header=TRUE)
aes<-read.csv("aes_delta_hyst.csv",header=TRUE)

 ###is there a correlation between drought and hysteranthy
#2003-4 and 2005-6 in germany Ivits et al 2014

##this plots the stations to show that they are all pretty mcuh in germany
newmap <- getMap(resolution = "low")
plot(newmap,
     xlim = c(0, 20),
     ylim = c(40, 55),
     asp = 1)
points(aln$lon,aln$lat,col="red")
points(frax$lon,frax$lat,col="blue")
points(aes$lon,aes$lat,col="purple")
points(15.28300,45.2670,col="green")###gerride of this point
###pretty much all site are in germany except for 1 alnus in croatia
aln<-dplyr::filter(aln,lat!=45.2670)

###subset to the 4 drought year
aln.sub<-dplyr::filter(aln,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
aln.sub$drought<-ifelse(aln.sub$year %in% c(2003,2004,2005,2006),1,0)

Anova(lmer(offset~drought+(1|s_id),data=aln.sub),type=3)
Anova(lmer(flower~drought+(1|s_id),data=aln.sub),type=3)

ggplot(aln.sub,aes(as.factor(drought),offset))+geom_boxplot()+ggtitle("Alnus incana")
ggplot(aln.sub,aes(as.factor(drought),flower))+geom_boxplot()+ggtitle("Alnus incana") 

frax.sub<-dplyr::filter(frax,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
frax.sub$drought<-ifelse(frax.sub$year %in% c(2003,2004,2005,2006),1,0)

Anova(lmer(offset~drought+(1|s_id),data=frax.sub))
Anova(lmer(flower~drought+(1|s_id),data=frax.sub))

ggplot(frax.sub,aes(as.factor(drought),offset))+geom_boxplot()+ggtitle("Fraxinus excelsior")
ggplot(frax.sub,aes(as.factor(drought),flower))+geom_boxplot()+ggtitle("Fraxinus excelsior")

aes.sub<-dplyr::filter(aes,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
aes.sub$drought<-ifelse(aes.sub$year %in% c(2003,2004,2005,2006),1,0)

Anova(lmer(offset~drought+(1|s_id),data=aes.sub))
Anova(lmer(flower~drought+(1|s_id),data=aes.sub))
ggplot(aes.sub,aes(as.factor(drought),offset))+geom_boxplot()+ggtitle("Aesculus")
ggplot(aes.sub,aes(as.factor(drought),flower))+geom_boxplot()+ggtitle("Aesculus")

###
library("raster")
library("remote")
library(reshape2)

d<-frax
calc<-d %>% group_by(s_id) %>% summarise(ave.offset=mean(offset))
calc2<-d %>% group_by(s_id) %>% summarise(ave.flower=mean(flower))
d<-dplyr::select(d,s_id,lon,lat)

d<-left_join(d,calc)
d<-left_join(d,calc2)
d<-d[!duplicated(d), ]

moist<-raster("grids_germany_multi_annual_soil_moist_1991-2010_05.asc") ##Gause Kruger 3


coordinates(d)<- ~lon + lat
proj4string(d) <- CRS("+proj=longlat +datum=WGS84")

p <- spTransform(d, CRS("+proj=tmerc +lat_0=50.625 +lon_0=9.84375, +k=1 +x_0=3559832.734474 +y_0=5610860.187573 +ellps=krass +units=m +no_defs"))
gaas<-coordinates(p)

goo<-as.data.frame(gaas)
goo2<-as.data.frame(d)
colnames(goo)<-c("x","y")

class(d)
goot<-cbind(goo2,goo)


###decent projectiontry to extract
foo<-extract(moist, matrix(c(goot$x,goot$y), ncol = 2))
goot$SM<-foo
class(goot)

summary(lm(ave.offset~SM,data=goot))


plot(moist)
points(p)

###alternative method find the average SM at the weather station closed to my points
library("rdwd")
nearbyStations(54.5000, 11.1,radius=5)
?nearbyStations()
moist

link <- selectDWD("Fehmarn", res="mutli_annual", var="soil", per="recent")
