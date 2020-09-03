rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/investment")
library("housingData")
library(stringr)
  sp<-read.csv("herbaria_prunus_rec.csv")


head(geoCounty)
head(sp$county)
geoCounty$rMapState<-str_to_title(geoCounty$rMapState)
colnames(geoCounty)[6]<-"stateProvince"

goo<-dplyr::left_join(sp,geoCounty,by=c("county","stateProvince"))

library("ncdf4")
library(raster)

nc_data<-nc_open("lbda-v2_kddm_pmdi_2017.nc")
print(nc_data)

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

pdsi.array <- ncvar_get(nc_data, "PDSI")

dim(pdsi.array)
fillvalue <- ncatt_get(nc_data, "PDSI", "_FillValue")
fillvalue
pdsi.array[pdsi.array == fillvalue$value] <- NA
pdsi.slice <- pdsi.array[, , 1] 
r <- raster(t(pdsi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))            
plot(pdsi.array)

palmer.b <- brick("lbda-v2_kddm_pmdi_2017.nc")  

lonpoints<-goo$lon
latpoints<-goo$lat
extract.pts <- cbind(lonpoints,latpoints)
ext <- extract(palmer.b,extract.pts,method="simple")

??extract()
plot(palmer.b)
points(lonpoints,latpoints,pch=4,col="red")

mean.goo <- calc(palmer.b, fun = mean,na.rm=TRUE)
ext <- extract(mean.goo,extract.pts,method="simple")

goo$pdsi<-ext

hyster<-read.csv("..//Data/rosaceae.csv")
head(hyster)
hyster<-dplyr::filter(hyster,genus=="Prunus")
colnames(hyster)[3]<-"specificEpithet"
goo2<-dplyr::left_join(goo,hyster,by="specificEpithet")
unique(goo$specificEpithet)
unique(hyster$specificEpithet)

goo2<-dplyr::filter(goo2,!is.na(hysteranthy))

png("pdsiboxes.png",width = 10,height=10,units = "in",res = 300)
ggplot(goo2,aes(as.factor(hysteranthy),pdsi))+geom_boxplot(outlier.shape = NA)+theme_minimal()+coord_cartesian(ylim = c(-.6, .6))
dev.off()
range(goo2$pdsi,na.rm = TRUE)
?stat_summary
car::Anova(lm(pdsi~as.factor(hysteranthy),data=goo2),type="III")

?geom_boxplot()
png("pdsimaps.png",width = 10,height=10,units = "in",res = 300)
plot(mean.goo)
points(lonpoints,latpoints, col=c("royal blue","black"),pch=c(8,5),cex=0.7)
legend(-80,20,c("hysteranthous","seranthous"),pch=c(8,5),col=c("royal blue","black"))
printa
dev.off()
