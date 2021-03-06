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
frax50<-read.csv("input/frax50year.csv")
aesc20<-read.csv("input/aesc20year.csv")
alnu20<-read.csv("input/alnu20year.csv")
alnu50<-read.csv("input/alnu50year.csv")
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
intra.df$soil.cent.neg<--(intra.df$soil.cent)

intra.df<-filter(intra.df, year>1991)
table(intra.df.hyst$taxa)

intra.df.hyst<-filter(intra.df, FLS>=0)
goo<-lm(FLS~soil.cent*flo.cent.neg,data=intra.df.hyst)
summary(goo)
rsqr(goo)

intra.df.ser<-filter(intra.df, FLS<=0)
table(intra.df.ser$taxa)
goo2<-lmer(FLS~soil.cent*flo.cent+(soil.cent*flo.cent|taxa),data=intra.df.ser)
summary(goo2)


df.intra.alnus<-filter(intra.df,taxa=="ALN.GLU")
df.intra.frax<-filter(intra.df,taxa=="FRA.EXC")
df.intra.bet<-filter(intra.df,taxa=="BET.PEN")
df.intra.aes<-filter(intra.df,taxa=="AES.HIP")
df.intra.fag<-filter(intra.df,taxa=="FAG.SYL")
df.intra.tili<-filter(intra.df,taxa=="TIL.HET")


frax.ind<-lmer(FLS~soil.cent+(soil.cent|year),data=df.intra.alnus)
summary(frax.ind)
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
library("ggrepel")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

alnumans<- frax50 %>% group_by(s_id,lat,lon) %>% summarise(meanFLS=mean(FLS),sdFLS=sd(FLS)) 
alnumans$meanFLS<-(floor(alnumans$meanFLS))
alnumans$sdFLS<-(floor(alnumans$sdFLS))
alnumans$hyst<-ifelse(alnumans$meanFLS>0,"hyst","ser")
alnumans$sdFLS <- paste0("(",alnumans$sdFLS , ")")
alnumans$FLS <- paste0(alnumans$meanFLS,alnumans$sdFLS)
randomRows <- function(df,n){
  return(df[sample(nrow(df),n),])
}
shortal<-randomRows(alnumans,20)


ggplot(data = world) + geom_sf()+
  geom_point(data = shortal, aes(x = lon, y = lat, fill=hyst), size = 1, 
             shape = 23) +
  coord_sf(xlim = c(5.5, 16), ylim = c(46, 55.5), expand = FALSE)

setEPS()
postscript("fraxmaps.eps",width = 7, height = 8)
ggplot(data = world) + geom_sf()+geom_text_repel(data = shortal, aes(x = lon, y = lat, label = FLS),size=4, 
                   nudge_x = c(.5, -.5, 1, 2, -1), nudge_y = c(0.25, -0.25, 0.5, 0.5, -0.5)) +
  geom_point(data = alnumans, aes(x = lon, y = lat), size = .5)+
  coord_sf(xlim = c(5.5, 16), ylim = c(47, 55), expand = TRUE)
dev.off()
frax.site<-lmer(FLS~soil.cent+flo.cent+(1|s_id),data=df.intra.frax)
summary(frax.site)

aln.f<-lm(FLS~flo.cent,data=df.intra.alnus)
aln.l<-lm(FLS~leaf.cent,data=df.intra.alnus)
frax.f<-lm(FLS~flo.cent,data=df.intra.frax)
frax.l<-lm(FLS~leaf.cent,data=df.intra.frax)
aes.f<-lm(FLS~flo.cent,data=df.intra.aes)
aes.l<-lm(FLS~leaf.cent,data=df.intra.aes)
bet.f<-lm(FLS~flo.cent,data=df.intra.bet)
bet.l<-lm(FLS~leaf.cent,data=df.intra.bet)
fag.f<-lm(FLS~flo.cent+leaf.cent,data=df.intra.fag)
fag.l<-lm(FLS~leaf.cent,data=df.intra.fag)
tili.f<-lm(FLS~flo.cent+leaf.cent,data=df.intra.tili)
tili.l<-lm(FLS~leaf.cent,data=df.intra.tili)

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


