rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/investment")
library("housingData")
library(stringr)
library("ncdf4")
library(raster)
library(ggplot2)
library("brms")

sp<-read.csv("herbaria_prunus_rec.csv")

allpru<-sp %>%group_by(specificEpithet)%>% count() %>% arrange(desc(n))

geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" ## make column names compatible
prunus.data<-dplyr::left_join(sp,geoCounty,by=c("county","stateProvince"))

geo.sp<-filter(prunus.data,county!="")
allpru.geo<-geo.sp %>%group_by(specificEpithet)%>% count() %>% arrange(desc(n))
all.specs<-merge(allpru,allpru.geo,by="specificEpithet")
colnames(all.specs)[c(2,3)]<-c("n","n.wcounty")

all.specs<-all.specs %>% arrange(desc(n.wcounty))
 ## This assigns each country coordinates



pruno<-c("alleghaniensis","angustifolia","americana" ,"gracilis","geniculata","hortulana" ,"maritima",
  "mexicana","murrayana","munsoniana","nigra","rivularis","umbellata","subcordata","texana" )


xtable::xtable(all.specs)

prunocerasus<-dplyr::filter(all.specs, specificEpithet %in% pruno)



palmer.b <- brick("lbda-v2_kddm_pmdi_2017.nc")  ## read in balmer drought index data 

plot(palmer.b[[1900:2000]])
lonpoints<-prunocerasus$lon # make vector of prunus coordinates
latpoints<-prunocerasus$lat #
extract.pts <- cbind(lonpoints,latpoints) #make coordinates to extract
ext <- extract(palmer.b,extract.pts,method="simple") ###extract drought info from each coordinate of prunus colllect



mean.prunus <- calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
ext <- extract(mean.prunus,extract.pts,method="simple")

prunocerasus$pdsi<-ext

par(mar=c(1,1,1,1))
plot(mean.prunus)
points(lonpoints,latpoints,cex=0.7)

hyster<-read.csv("cherry_data.csv") ## read in flower size, phenology and fruiot size data


colnames(hyster)[1]<-"specificEpithet"
hyster$mean_size<-(hyster$petal_high+hyster$petal_low)/2
hyster$mean_num<-(hyster$inflor_high+hyster$inflor_low)/2
hyster$mean_fruit<-(hyster$fruit_high+hyster$fruit_low)/2
hyster$mean_axis<-(hyster$midrib_high+hyster$midrib_low)/2
hyster$mean_ped<-(hyster$pedicel_high+hyster$pedicel_low)/2
hyster$display<-hyster$mean_num*hyster$mean_size
hyster$hyst<-ifelse(hyster$FLS %in%c("before","before/with"),"hys","ser")
hyster$hyst2<-ifelse(hyster$FLS %in%c("before"),"hys","ser")
head(hyster)
pruno<-c("alleghaniensis","angustifolia","americana" ,"gracilis","geniculata","hortulana" ,"maritima",
         "mexicana","murrayana","munsoniana","nigra","rivularis","umbellata","subcordata","texana" )

prunycer<-dplyr::filter(hyster,specificEpithet %in% pruno)
prunycer2<-left_join(prunycer,prunocerasus)
#d<-dplyr::filter(prunus.data.2,!is.na(FLS)) ### 7668 rows have county coordiates

ggplot(prunycer,aes(FLS,mean_size))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")

ggplot(prunycer,aes(hyst,mean_size))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")

ggplot(prunycer,aes(hyst2,mean_size))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")

ggplot(prunycer2,aes(FLS,pdsi))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")


