rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

setwd("~/Documents/git/proterant/sub_projs/PEP725/")
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(tibble)
library(raster)
library(lme4)

frax20<-read.csv("input/frax20year.csv")
aesc20<-read.csv("input/aesc20year.csv")
alnu20<-read.csv("input/alnu20year.csv")
fagu20<-read.csv("input/fagu20year.csv")
tili20<-read.csv("input/tili20year.csv")
betu20<-read.csv("input/betu20year.csv")

frax20$taxa<-"FRA.EXC"
aesc20$taxa<-"AES.HIP"
alnu20$taxa<-"ALN.GLU"
fagu20$taxa<-"FAG.SYL"
tili20$taxa<-"TIL.HET"
betu20$taxa<-"BET.PEN"

fullcomps<-rbind(frax20,aesc20,alnu20,fagu20,tili20,betu20)

moist.aug<-raster("climate_data/grids_germany_multi_annual_soil_moist_1991-2010_08.asc") ##Gause Kruger 3 for August
moist.apr<-raster("climate_data/grids_germany_multi_annual_soil_moist_1991-2010_04.asc")

coordinates(fullcomps)<- ~lon + lat ###PEP site coordinates
proj4string(fullcomps) <- CRS("+proj=longlat +datum=WGS84") ##establish the project


p <- spTransform(fullcomps, CRS("+proj=tmerc +lat_0=50.625 +lon_0=9.84375, +k=1 +x_0=3559832.734474 +y_0=5610860.187573 +ellps=krass +units=m +no_defs")) ###convert our coordinate to Gausse Kruger
gaas<-coordinates(p) ### make these values coordinates in structure

gaas.cord<-as.data.frame(gaas)
pep.cord<-as.data.frame(fullcomps)
colnames(gaas.cord)<-c("x","y")
intra.df<-cbind(pep.cord,gaas.cord) 

Soil<-extract(moist.aug, matrix(c(intra.df$x,intra.df$y), ncol = 2)) ### extract the soil moisture at the pep sites
intra.df$SM<-Soil

intra.df$flo.cent<-(intra.df$flo_day-mean(intra.df$flo_day,na.rm=TRUE))/(sd(intra.df$flo_day,na.rm=TRUE))
intra.df$flo.cent.neg<--(intra.df$flo_day-mean(intra.df$flo_day,na.rm=TRUE))/(sd(intra.df$flo_day,na.rm=TRUE))
intra.df$leaf.cent<-(intra.df$leaf_day-mean(intra.df$leaf_day,na.rm=TRUE))/(sd(intra.df$leaf_day,na.rm=TRUE))
intra.df$soil.cent<-(intra.df$SM-mean(intra.df$SM,na.rm=TRUE))/(sd(intra.df$SM,na.rm=TRUE))

wind<-c("ALN.GLU")
df.intra.alnus<-filter(intra.df,taxa=="ALN.GLU")
df.intra.frax<-filter(intra.df,taxa=="FRA.EXC")
df.intra.bet<-filter(intra.df,taxa=="BET.PEN")
df.intra.aes<-filter(intra.df,taxa=="AES.HIP")
df.intra.fag<-filter(intra.df,taxa=="FAG.SYL")
df.intra.tili<-filter(intra.df,taxa=="TIL.HET")


###does fleaf time or flow time predict fls

####soil moisture by species
aln.mod<-lm(FLS~soil.cent*flo.cent.neg,data=df.intra.alnus)
summary(aln.mod)

bet.mod<-lm(FLS~soil.cent*flo.cent.neg,data=df.intra.bet)
summary(bet.mod)

aes.mod<-lm(FLS~soil.cent*flo.cent.neg,data=df.intra.aes)
summary(aes.mod)

tili.mod<-lm(FLS~soil.cent*flo.cent.neg,data=df.intra.tili)
summary(tili.mod)

frax.mod<-lm(FLS~soil.cent*flo.cent.neg,data=df.intra.frax)
summary(frax.mod)
##3


