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

#read in clean data of offset for each species
aln<-read.csv("alnus10_delta_hyst.csv",header=TRUE)
frax<-read.csv("fraxinus10_delta_hyst.csv",header=TRUE)
#aes<-read.csv("aes_delta_hyst.csv",header=TRUE)
quer<-read.csv("Quercus10_delta_hyst.csv",header=TRUE)
pru<-read.csv("Prunus10_delta_hyst.csv",header=TRUE)
cory<-read.csv("cory10_delta_hyst.csv",header=TRUE) 
bet<-read.csv("betpen10_delta_hyst.csv",head=TRUE)
###is there a correlation between drought and hysteranthy
#2003-4 and 2005-6 in germany Ivits et al 2014

##this plots the stations to show that they are all pretty mcuh in germany
#newmap <- getMap(resolution = "low")
#plot(newmap,
 #    xlim = c(9, 11),
  #   ylim = c(47, 55),
   #  asp = 1)
#points(aln$lon,aln$lat,col="red")
#points(frax$lon,frax$lat,col="blue")
#points(aes$lon,aes$lat,col="purple")
#points(15.28300,45.2670,col="green")###gerride of this point
###pretty much all site are in germany except for 1 alnus in croatia
aln<-dplyr::filter(aln,lat!=45.2670)

#######################################################################################
###PART I, Does Hysteranthy chance in a DROUGHT YEAR
###subset to the 4 drought year, 4 wet years
aln.sub<-dplyr::filter(aln,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
aln.sub$drought<-ifelse(aln.sub$year %in% c(2003,2004,2005,2006),1,0) #1 are drought

Anova(lmer(offset~drought+(1|s_id),data=aln.sub),type=3)
summary(lmer(offset~drought+(1|s_id),data=aln.sub),type=3)##drought reduces offset
summary(lmer(flower~drought+(1|s_id),data=aln.sub),type=3) ##drought delays flowerinf

ggplot(aln.sub,aes(as.factor(drought),offset))+geom_boxplot()+ggtitle("Alnus incana")
ggplot(aln.sub,aes(as.factor(drought),flower))+geom_boxplot()+ggtitle("Alnus incana") 

frax.sub<-dplyr::filter(frax,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
frax.sub$drought<-ifelse(frax.sub$year %in% c(2003,2004,2005,2006),1,0)

summary(lmer(offset~drought+(1|s_id),data=frax.sub)) ##no effect of drought on alnus
summary(lmer(flower~drought+(1|s_id),data=frax.sub)) ##drought delays flowerig

ggplot(frax.sub,aes(as.factor(drought),offset))+geom_boxplot()+ggtitle("Fraxinus excelsior")
ggplot(frax.sub,aes(as.factor(drought),flower))+geom_boxplot()+ggtitle("Fraxinus excelsior")

aes.sub<-dplyr::filter(aes,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
aes.sub$drought<-ifelse(aes.sub$year %in% c(2003,2004,2005,2006),1,0)

summary(lmer(offset~drought+(1|s_id),data=aes.sub)) ##no effect on offset
summary(lmer(flower~drought+(1|s_id),data=aes.sub)) ##drought delays flowering
ggplot(aes.sub,aes(as.factor(drought),offset))+geom_boxplot()+ggtitle("Aesculus")
ggplot(aes.sub,aes(as.factor(drought),flower))+geom_boxplot()+ggtitle("Aesculus")

#################################NO################################################################

###################################################################################
#PART II Does soil average moisture predict hysteranthy?###########################
library("raster")
library("remote")
library(reshape2)
aln$taxa<-"A. glutinosa"
frax$taxa<-"F. excelsior"
  
d<-rbind(aln,frax,cory,bet,quer,pru) #make a datasheet of fraxinus and alnus
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


###decent projectiontry to extract
Soil<-extract(moist, matrix(c(goot$x,goot$y), ncol = 2))
goot$SM<-Soil
class(goot)

droughtmod<-brm(ave.offset~SM+taxa+SM:taxa,data=goot)
nosp<-brm(ave.offset~SM,data=goot)
#coef(droughtmod)
summary(droughtmod) ### as SM decreases offset increases (if I did this right)
summary(nosp)

SMresults<-as.data.frame(tidy(droughtmod,robust = TRUE))
SMresults<-SMresults %>% "["(.,1:6,)
SMresults$term[SMresults$term=="b_SM"]<-"Soil moisture"
SMresults$term[SMresults$term=="b_Intercept"]<-"Intercept"
SMresults$term[SMresults$term=="b_taxaBetulapendula"]<-"B. Pendula"
SMresults$term[SMresults$term=="b_taxaF.excelsior"]<-"F. excelsior"
SMresults$term[SMresults$term=="b_SM:taxaF.excelsior"]<-"soil:Franxinus"
SMresults$term[SMresults$term=="b_SM:taxaBetulapendula"]<-"soil:Betula"

unique(SMresults$term)
zoom<-SMresults %>% "["(.,5:6,)
zoom1<-SMresults %>% "["(.,2:2,)
zoom<-rbind(zoom,zoom1)
pd=position_dodgev(height=0.4)
X<-ggplot(SMresults,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,,width=.1)+geom_vline(aes(xintercept=0),color="black")+xlim(-50,50)
Y<-ggplot(zoom,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,width=.1)+geom_vline(aes(xintercept=0),color="black")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())+xlim(-0.2,0.2)

vp <- viewport(width = 0.46, height = 0.53, x = .32, y = 1,just=c("left","top"))
X
print(Y, vp = vp )  


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
SMFresults$term[SMresults$term=="b_Intercept"]<-"Intercept"
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


vp <- viewport(width = 0.46, height = 0.53, x = .32, y = 1,just=c("left","top"))
X
print(Y, vp = vp )  


flo<-lm(ave.flower~ave_last+Taxa+Taxa:ave_last,data=last.freeze.data)###last freeze doesn't matter either
leaf<-lm(ave.leaf~ave_last+Taxa+Taxa:ave_last,data=last.freeze.data)###last freeze doesn't matter either
##last freeze doesn't matter either
summary(freezmod)
summary(flo)
summary(leaf)

freeze.results<-as.data.frame(tidy(freezmod,robust = TRUE))
freeze.results<-freeze.results %>% "["(.,1:3,)
zoom<-freeze.results<-freeze.results %>% "["(.,2:3,)
pd=position_dodgev(height=0.4)
X<-ggplot(freeze.results,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,,width=.1)+geom_vline(aes(xintercept=0),color="black")+theme_base()
Y<-ggplot(zoom,aes(estimate,term))+geom_point(position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper),position=pd,width=.1)+geom_vline(aes(xintercept=0),color="black")+theme(axis.title.y = element_blank(),axis.text.y = element_blank())

vp <- viewport(width = 0.5, height = 0.5, x = 0.45, y = .65,just=c("left","top"))
X
print(Y, vp = vp, )                                               
                                           
?viewport()
