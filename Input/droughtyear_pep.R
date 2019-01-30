rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(tidyverse)
library(rworldmap)
library(lme4)
library(car)
library(rgdal)
library(ggstance)
library(grid)
library(ggthemes)
library(brms)
library(lubridate)
library(broom)
library("ggmap")
library(sp)
library("raster")
library("remote")
library(reshape2)


#load("intervar.RData")

#read in clean data of offset for each species
aln<-read.csv("alnus10_delta_hyst.csv",header=TRUE)
aln2<-read.csv("alnus_delta_hyst.csv")
frax<-read.csv("fraxinus10_delta_hyst.csv",header=TRUE)
frax2<-read.csv("fraxinus_delta_hyst.csv")
aes<-read.csv("aes10_delta_hyst.csv",header=TRUE)
aes2<-read.csv("aes_delta_hyst.csv",header=TRUE)
quer<-read.csv("Quercus10_delta_hyst.csv",header=TRUE)
pru<-read.csv("Prunus10_delta_hyst.csv",header=TRUE)
cory<-read.csv("cory10_delta_hyst.csv",header=TRUE) 
bet<-read.csv("betpen10_delta_hyst.csv",head=TRUE)


aln$taxa<-"Alnus glutinosa"
frax$taxa<-"Fraxinus excelsior"
aes$taxa<-"Aesculus hippocastenum"

d<-rbind(aln,frax,cory,bet,quer) #make a datasheet of fraxinus and alnus
d<-filter(d,!is.na(offset))

###read in drought indexmap

DI.aug<-raster("grids_germany_multi_annual_drought_index_1981-2010_08.asc") 
plot(DI.aug)
DI.apr<-raster("grids_germany_multi_annual_drought_index_1981-2010_04.asc")
plot(DI.apr)
DI.sept<-raster("grids_germany_multi_annual_drought_index_1981-2010_09.asc")
plot(DI.sept)

moist<-raster("grids_germany_multi_annual_soil_moist_1991-2010_08.asc") ##Gause Kruger 3 for August
plot(moist)

moist2<-raster("grids_germany_multi_annual_soil_moist_1991-2010_04.asc")


###auge seems to have the most variation

coordinates(d)<- ~lon + lat ###PEP site coordinates
proj4string(d) <- CRS("+proj=longlat +datum=WGS84")

# THis should be the rihgt but doesnt project p <- spTransform(d, CRS("+proj=tmerc +lat_0=47.2700 +lon_0=7.5000 +k=1 +x_0=3386564.9400+y_0=5237917.9109 +ellps=krass +units=m +no_defs"))
p <- spTransform(d, CRS("+proj=tmerc +lat_0=50.625 +lon_0=9.84375, +k=1 +x_0=3559832.734474 +y_0=5610860.187573 +ellps=krass +units=m +no_defs"))

gaas<-coordinates(p) ### p transforms site to weird german gauss projection

gaas.cord<-as.data.frame(gaas)
pep.cord<-as.data.frame(d)
colnames(gaas.cord)<-c("x","y")

class(d)
df<-cbind(pep.cord,gaas.cord)


###decent projectiontry to extract
Soil<-extract(moist, matrix(c(df$x,df$y), ncol = 2))
df$SM<-Soil

drought.index<-extract(DI.aug, matrix(c(DI.df$x,DI.df$y), ncol = 2))

df$drought_index<-drought.index
###
plot(moist)
points(df$x,DI.df$y,pch=16)


DI.mod<-lmer(offset~SM+(1|taxa),data=df) ###sign flip when you add flowering
summary(DI.mod)
coef(DI.mod)
AIC(DI.mod)

FLO.MOD<-lmer(offset~flower+(1|taxa),data=df)
summary(FLO.MOD)
coef(FLO.MOD)

flo.sm<-lmer(offset~flower+SM+(1|taxa,data=df))
summary(flo.sm)

                          
             ## #is there a correlation between drought and hysteranthy
#2003-4 and 2005-6 in germany Ivits et al 2014





amean.offset<- aln2 %>% group_by(s_id,lat,lon) %>% summarise(mean.offset=mean(offset))
amean.offset$mean.offset<-ceiling(amean.offset$mean.offset)

sd.offset<- aln2 %>% group_by(s_id,lat,lon) %>% summarise(sd.offset=sd(offset))
sd.offset$sd.offset<-ceiling(sd.offset$sd.offset)

stat<-left_join(amean.offset,sd.offset)
stat$value<-NA
stat$value <- paste(stat$mean.offset, stat$sd.offset, sep=",")


newmap <- getMap(resolution = "low")
plot(newmap,
     xlim = c(9.5, 10),
     ylim = c(47.5, 55),
     asp = 1)
goo<-aov(aln$offset~as.factor(aln$s_id))
TukeyHSD(goo)

points(stat$lon,stat$lat,col="blue",pch=16,cex=0.6)

text(stat$lon, y = stat$lat, stat$value, pos = 3, cex=0.6)



#points(frax$lon,frax$lat,col="blue")
#points(aes$lon,aes$lat,col="purple")
#points(15.28300,45.2670,col="green")###gerride of this point
###pretty much all site are in germany except for 1 alnus in croatia
aln<-dplyr::filter(aln,lat!=45.2670)

#######################################################################################
###PART I, Does Hysteranthy chance in a DROUGHT YEAR
###subset to the 4 drought year, 4 wet years
d.sub<-dplyr::filter(d,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
d.sub$drought<-ifelse(d.sub$year %in% c(2003,2004,2005,2006),1,0) #1 are drought


annual<-lmer(offset~drought*taxa+(1|s_id),dat=d.sub)
summary(annual)


#################################NO################################################################

###################################################################################
#PART II Does soil average moisture predict hysteranthy?###########################
library("raster")
library("remote")
library(reshape2)
aln$taxa<-"Alnus glutinosa"
frax$taxa<-"Fraxinus excelsior"
aes$taxa<-"Aesculus hippocastenum"
  
d<-rbind(aln,frax,cory,bet,quer,pru,aes) #make a datasheet of fraxinus and alnus
d<-filter(d,!is.na(offset))




calc<-d %>% group_by(s_id,taxa) %>% summarise(ave.offset=mean(offset)) ##average offset/sp/site
calc2<-d %>% group_by(s_id,taxa) %>% summarise(ave.flower=mean(flower)) ##average flowertime/sp/site
calc3<- d %>% group_by(s_id,taxa) %>% summarise(ave.leaf=mean(leaf)) ##average leaftime /sp/site
d<-dplyr::select(d,s_id,lon,lat,taxa) ### reudce datashet

d<-left_join(d,calc)
d<-left_join(d,calc2)
d<-left_join(d,calc3)
d<-d[!duplicated(d), ] ###join all datasheet so each site has an a3 average phoenogy parenpers for each sp
table(d$taxa)

newmap <- getMap(resolution = "low")
plot(newmap,
  xlim = c(9, 11),
  ylim = c(47, 60),
 asp = 1)
points(d$lon,d$lat,col="red")

moist<-raster("grids_germany_multi_annual_soil_moist_1991-2010_08.asc") ##Gause Kruger 3 for August
moist2<-raster("grids_germany_multi_annual_soil_moist_1991-2010_04.asc")
plot(moist)
#plot(moist2)
coordinates(d)<- ~lon + lat ###PEP site coordinates
proj4string(d) <- CRS("+proj=longlat +datum=WGS84") #and theith projection

p <- spTransform(d, CRS("+proj=tmerc +lat_0=50.625 +lon_0=9.84375, +k=1 +x_0=3559832.734474 +y_0=5610860.187573 +ellps=krass +units=m +no_defs"))
gaas<-coordinates(p) ### p transforms site to weird german gauss projection

gaas.cord<-as.data.frame(gaas)
pep.cord<-as.data.frame(d)
colnames(gaas.cord)<-c("x","y")

class(d)
goot<-cbind(pep.cord,gaas.cord)




droughtmod<-brm(ave.offset~SM+ave.flower+taxa+SM:taxa+ave.flower:taxa,data=goot)
summary(droughtmod)

nosp<-brm(ave.offset~SM*ave.flower,data=goot)
#coef(droughtmod)
summary(nosp) ### as SM decreases offset increases (if I did this right)

goot$sm.z<-(goot$SM-mean(goot$SM,na.rm=TRUE))/sd(goot$SM,na.rm=TRUE)
goot$flo.z<-(goot$ave.flower-mean(goot$ave.flower,na.rm=TRUE))/sd(goot$ave.flower,na.rm=TRUE)

### each species:
fraxy<-filter(goot, taxa=="Fraxinus excelsior")
alny<-filter(goot, taxa=="Alnus glutinosa")
betty<-filter(goot, taxa=="Betula pendula")
aesy<-filter(goot, taxa=="Aesculus hippocastenum")

frax.sm.only<-brm(ave.offset~SM,data=fraxy)
droughtfrax<-brm(ave.offset~SM*ave.flower,data=fraxy)
droughtfrax.z<-brm(ave.offset~sm.z*flo.z,data=fraxy)

droughtaln<-brm(ave.offset~SM*ave.flower,data=alny)
droughtbet<-brm(ave.offset~SM*ave.flower,data=betty)
droughtbet.noint<-brm(ave.offset~SM+ave.flower,data=betty)
droughtaesy<-brm(ave.offset~SM*ave.flower,data=aesy)
droughtaesy.z<-brm(ave.offset~sm.z*flo.z,data=aesy)
droughtfrax.z<-brm(ave.offset~sm.z*flo.z,data=fraxy)
droughtaln.z<-brm(ave.offset~sm.z*flo.z,data=alny)
droughtbet.z<-brm(ave.offset~sm.z*flo.z,data=betty)
summary(frax.sm.only)
summary(droughtfrax.z)
summary(droughtaln.z)
summary(droughtaesy.z)
summary(droughtbet.z)
summary(droughtaln)
summary(droughtfrax)
summary(droughtbet)
summary(droughtbet.noint)
table(goot$taxa)



SMresults<-as.data.frame(tidy(droughtmod,robust = TRUE))
SMresults<-SMresults %>% "["(.,1:9,)
unique(SMresults$term)
SMresults$term[SMresults$term=="b_SM"]<-"Soil moisture"
SMresults$term[SMresults$term=="b_Intercept"]<-"a_Intercept"
SMresults$term[SMresults$term=="b_ave.flower"]<-"flowering time"
SMresults$term[SMresults$term=="b_taxaBetulapendula"]<-"B. Pendula"
SMresults$term[SMresults$term=="b_taxaFraxinusexcelsior"]<-"F. excelsior"
SMresults$term[SMresults$term=="b_SM:taxaFraxinusexcelsior"]<-"soil:Franxinus"
SMresults$term[SMresults$term=="b_SM:taxaBetulapendula"]<-"soil:Betula"
SMresults$term[SMresults$term=="b_ave.flower:taxaFraxinusexcelsior"]<-"flowering:Franxinus"
SMresults$term[SMresults$term=="b_ave.flower:taxaBetulapendula"]<-"flowering:Betula"

unique(SMresults$term)
zoom<-SMresults %>% "["(.,6:9,)
zoom1<-SMresults %>% "["(.,2:3,)
zoom<-rbind(zoom,zoom1)
pd=position_dodgev(height=0.4)
X<-ggplot(SMresults,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,,width=.01)+geom_vline(aes(xintercept=0),color="black")+xlim(-50,80)
Y<-ggplot(zoom,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,width=.01)+geom_vline(aes(xintercept=0),color="black")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())+xlim(-0.8,0.5)


vp <- viewport(width = 0.46, height = 0.67, x = .2, y = 1,just=c("left","top"))
jpeg("..//figure/SMxflo.jpeg")
X
print(Y, vp = vp )  
dev.off()



flo.mod<-brm(ave.offset~ave.flower+taxa+ave.flower:taxa,data=goot)
summary(flo.mod)


plot(moist)
points(gaas.cord)

###summary for neither alnus or fraxinus to we see an association between increased hysteranthy and drought

r<-brick("~/Desktop/tn_0.25deg_reg_v17.0.nc", varname="tn", sep="")


values <- raster::extract(r,d)

dclim <- cbind.data.frame(coordinates(d),values)
colnames(dclim)
dx<-melt(dclim, id.vars=c("lon","lat"))
head(dx)
dx<-dx%>%
  rename(date=variable)%>%
  rename(Tlowavg=value)

dx$date<-substr(dx$date, 2,11)
dx$Date<- gsub("[.]", "-", dx$date)
####


dx$month<-substr(dx$Date, 6, 7)
dx$month<-as.numeric(dx$month)
head(dx)
dx$year<-substr(dx$Date, 1,4)
dx$year<-as.numeric(dx$year)

dx$winter<-ifelse(dx$month>=9 | dx$month<=4, "winter", 0)
winter<-dx[(dx$winter=="winter"),]
winter<-winter[!is.na(winter$Tlowavg),]

dx$spring<-ifelse(dx$month>=2 & dx$month<=4, "spring", 0)
ddx<-dx[(dx$spring=="spring"),]
ddx<-ddx[!is.na(ddx$Tlowavg),]

ddx$year<-as.numeric(substr(ddx$Date, 0, 4))

ddx$lat.long<-paste(ddx$lat, ddx$lon)
#ddx$mat<-ave(ddx$Tavg, ddx$year, ddx$lat.long)
head(ddx)
####frrezing days

freezy<-filter(ddx,Tlowavg<=-2.2)
unique(freezy$month)
freezy$DOY<-yday(freezy$Date)
head(freezy)
lastfreeze<-freezy %>% group_by(lat.long,year) %>% summarise(last_freez=last(DOY))

head(lastfreeze)
lastfreeze<-lastfreeze %>% group_by(lat.long) %>% summarise(ave_last=mean(last_freez))
range(lastfreeze$ave_last)

goot$lat.long<-paste(goot$lat, goot$lon)
head(goot)
colnames(goot)
colnames(lastfreeze)
last.freeze.data<-left_join(goot,lastfreeze,by="lat.long")
colnames(last.freeze.data)
unique(last.freeze.data$Taxa)
alnus.only<-filter(last.freeze.data,Taxa=="A. glutinosa")
Frax.only<-filter(last.freeze.data,Taxa=="F. excelsior")

freezmod<-lm(ave.offset~ave_last+taxa+taxa:ave_last,data=last.freeze.data)###last freeze doesn't matter either
summary(freezmod)

freezeySM<-brm(ave.offset~ave_last+SM+taxa+taxa:ave_last+taxa:SM,data=last.freeze.data)

summary(freezeySM)
SMFresults<-as.data.frame(tidy(freezeySM,robust = TRUE))
SMFresults<-SMFresults %>% "["(.,1:9,)
SMFresults$term[SMFresults$term=="b_SM"]<-"soil moisture"
SMFresults$term[SMFresults$term=="b_ave_last"]<-"last frost"
SMFresults$term[SMFresults$term=="b_Intercept"]<-"B_Intercept"
SMFresults$term[SMFresults$term=="b_taxaBetulapendula"]<-"B. Pendula"
SMFresults$term[SMFresults$term=="b_taxaF.excelsior"]<-"F. excelsior"
SMFresults$term[SMFresults$term=="b_SM:taxaF.excelsior"]<-"soil:Franxinus"
SMFresults$term[SMFresults$term=="b_SM:taxaBetulapendula"]<-"soil:Betula"
SMFresults$term[SMFresults$term=="b_ave_last:taxaBetulapendula"]<-"frost:Betula"
SMFresults$term[SMFresults$term=="b_ave_last:taxaF.excelsior"]<-"frost:Fraxinus"
unique(SMFresults$term)

zoom<-SMFresults %>% "["(.,2:3,)
zoom1<-SMFresults %>% "["(.,6:9,)
zoom<-rbind(zoom,zoom1)
pd=position_dodgev(height=0.4)
X<-ggplot(SMFresults,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,,width=.1)+geom_vline(aes(xintercept=0),color="black")+xlim(-55,55)
Y<-ggplot(zoom,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,width=.1)+geom_vline(aes(xintercept=0),color="black")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+xlim(-0.2,0.2)


vp <- viewport(width = 0.46, height = 0.7, x = .32, y = 1,just=c("left","top"))
X
print(Y, vp = vp )  

save.image("intervar.RData")
