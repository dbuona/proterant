rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(tidyverse)
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
aln<-read.csv("alnus10_delta_hyst.csv",header=TRUE) ## alnus from stations with atleast 10 years of data
aln2<-read.csv("alnus_delta_hyst.csv") ## # with 50 years
frax<-read.csv("fraxinus10_delta_hyst.csv",header=TRUE)## 10
frax2<-read.csv("fraxinus_delta_hyst.csv") ## at least 50
aes<-read.csv("aes10_delta_hyst.csv",header=TRUE) ##10 (Non hysteramthous)
aes2<-read.csv("aes_delta_hyst.csv",header=TRUE) ## 50 years
quer<-read.csv("Quercus10_delta_hyst.csv",header=TRUE) ##10
pru<-read.csv("Prunus10_delta_hyst.csv",header=TRUE)##10
cory<-read.csv("cory10_delta_hyst.csv",header=TRUE)  ##10
bet<-read.csv("betpen10_delta_hyst.csv",head=TRUE) #10

unique

aln$taxa<-"Alnus glutinosa" ##assign species
frax$taxa<-"Fraxinus excelsior"
aes$taxa<-"Aesculus hippocastenum"

d<-rbind(aln,frax,cory,bet,quer) #make a datasheet of fraxinus and alnus
#d<-filter(d,!is.na(offset))

###read in drought indexmap

DI.aug<-raster("grids_germany_multi_annual_drought_index_1981-2010_08.asc") 
DI.apr<-raster("grids_germany_multi_annual_drought_index_1981-2010_04.asc")
DI.sept<-raster("grids_germany_multi_annual_drought_index_1981-2010_09.asc")
moist.aug<-raster("grids_germany_multi_annual_soil_moist_1991-2010_08.asc") ##Gause Kruger 3 for August
moist.apr<-raster("grids_germany_multi_annual_soil_moist_1991-2010_04.asc")
##use august ##Question 1: is this reasonable or should I do some sort of model selection?

coordinates(d)<- ~lon + lat ###PEP site coordinates
proj4string(d) <- CRS("+proj=longlat +datum=WGS84") ##establish the project

#   <- spTransform(d, CRS("+proj=tmerc +lat_0=47.2700 +lon_0=7.5000 +k=1 +x_0=3386564.9400+y_0=5237917.9109 +ellps=krass +units=m +no_defs")) ### this doesnt work
p <- spTransform(d, CRS("+proj=tmerc +lat_0=50.625 +lon_0=9.84375, +k=1 +x_0=3559832.734474 +y_0=5610860.187573 +ellps=krass +units=m +no_defs")) ###convert our coordinate to Gausse Kruger


gaas<-coordinates(p) ### make these values coordinates in structure

gaas.cord<-as.data.frame(gaas)
pep.cord<-as.data.frame(d)
colnames(gaas.cord)<-c("x","y")


df<-cbind(pep.cord,gaas.cord) ### aff the next


###decent projectiontry to extract
Soil<-extract(moist.aug, matrix(c(df$x,df$y), ncol = 2)) ### extract the soil moisture at the pep sites
df$SM<-Soil
drought.index<-extract(DI.aug, matrix(c(df$x,df$y), ncol = 2)) ### same but form drought index
df$drought_index<-drought.index
###
df$flo.cent<-df$flower-mean(df$flower,na.rm=TRUE)
df$leaf.cent<-df$leaf-mean(df$leaf,na.rm=TRUE)
df$soil.cent<-df$SM-mean(df$SM,na.rm=TRUE)


unique(df$taxa)
DI.mod<-lm(offset~SM*taxa,data=df) 
summary(DI.mod) ### seems like alnus is the only 

df.alnus<-filter(df,taxa=="Alnus glutinosa")
DI.mod.aln<-lm(offset~SM,data=df.alnus) 

df.frax<-filter(df,taxa=="Fraxinus excelsior")
DI.mod.frax<-lm(offset~SM,data=df.frax) 
summary(DI.mod.frax)

df.bet<-filter(df,taxa=="Betula pendula")
DI.mod.bet<-lm(offset~SM,data=df.bet) 
summary(DI.mod.bet) ####becomse more leafing first
####





leaf.aln<-lm(offset~leaf.cent,data=df.alnus)
summary(leaf.aln)
flo.aln<-lm(offset~flo.cent,data=df.alnus)
summary(flo.aln)
soil.aln<-lm(offset~soil.cent,data=df.alnus) 
summary(soil.aln)

twopred.aln<-lm(offset~flo.cent*soil.cent,data=df.alnus)
summary(twopred.aln)
twopred.frax<-lm(offset~flo.cent*soil.cent,data=df.frax)
summary(twopred.frax)
twopred.bet<-lm(offset~flo.cent*soil.cent,data=df.bet)
summary(twopred.bet)



###### Other question, are propualtions even "significantly" different
aov.aln<-lm(offset~as.factor(s_id), data=aln2)
Anova(aov.aln,ype="III")
###yest



#######################################################################################
###PART I, Does Hysteranthy chance in a DROUGHT YEAR
###subset to the 4 drought year, 4 wet years

dat<-rbind(aln,frax,cory,bet,quer)
d.sub<-dplyr::filter(dat,year %in% c(2003,2004,2005,2006,2007,2008,2009,2010))
d.sub$drought<-ifelse(d.sub$year %in% c(2003,2004,2005,2006),1,0) #1 are drought


annual<-lmer(offset~drought+(1|s_id)+(drought|taxa),dat=d.sub) ### drought years have less ofset
summary(annual)
coef(annual)

annual.flower<-lmer(flower~drought+(1|s_id)+(drought|taxa),dat=d.sub) ###droguht delayes flowering
coef(annual.flower)

annual.leaf<-lmer(leaf~drought+(1|s_id)+(drought|taxa),dat=d.sub) 
coef(annual.leaf)

####should I do this with soil annual soil moisture data at the sites instead
#################################NO################################################################

###################################################################################

#####do place with later frost affect offset ####probably not. 
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


stop()##################################################################
######plot
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



points(stat$lon,stat$lat,col="blue",pch=16,cex=0.6)

text(stat$lon, y = stat$lat, stat$value, pos = 3, cex=0.6)



#points(frax$lon,frax$lat,col="blue")
#points(aes$lon,aes$lat,col="purple")
#points(15.28300,45.2670,col="green")###gerride of this point
###pretty much all site are in germany except for 1 alnus in croatia



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

