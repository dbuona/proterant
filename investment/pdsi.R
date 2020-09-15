rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/investment")
library("housingData")
library(stringr)
library("ncdf4")
library(raster)

sp<-read.csv("herbaria_prunus_rec.csv") ## all the prunus records


geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" ## make column names compatible

prunus.data<-dplyr::left_join(sp,geoCounty,by=c("county","stateProvince")) ## This assigns each country coordinates


palmer.b <- brick("lbda-v2_kddm_pmdi_2017.nc")  ## read in balmer drought index data 

lonpoints<-prunus.data$lon # make vector of prunus coordinates
latpoints<-prunus.data$lat #
extract.pts <- cbind(lonpoints,latpoints) #make coordinates to extract
ext <- extract(palmer.b,extract.pts,method="simple") ###extract drought info from each coordinate of prunus colllect



mean.prunus <- calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
ext <- extract(mean.prunus,extract.pts,method="simple")

prunus.data$pdsi<-ext ##3 noiw your data has a pdsi estimate for each county of collect

hyster<-read.csv("..//Data/rosaceae.csv") ## read in flower size, phenology and fruiot size data

hyster<-dplyr::filter(hyster,genus=="Prunus") ## subset to jsut prunus
colnames(hyster)[3]<-"specificEpithet"
hyster$mean_size<-(hyster$flower_size_high+hyster$flos_per_low)/2
hyster$mean_num<-(hyster$flos_per_high+hyster$flos_per_low)/2
hyster$display<-hyster$mean_num*hyster$mean_size


prunus.data.2<-dplyr::left_join(prunus.data,hyster,by="specificEpithet")


d<-dplyr::filter(prunus.data.2,!is.na(hysteranthy)) ### 7668 rows have county coordiates




## hypothesis1: hysteranthy relates to flower size

ggplot(hyster,aes(mean_size))+geom_histogram(bins=20)+facet_wrap(~as.factor(hysteranthy))+geom_vline(aes(xintercept=mean(mean_size),color="red"))+geom_vline(aes(xintercept=median(mean_size),color="blue"))
ggplot(hyster,aes(mean_num))+geom_histogram(bins=20)+facet_wrap(~as.factor(hysteranthy))+geom_vline(aes(xintercept=mean(mean_num),color="red"))+geom_vline(aes(xintercept=median(mean_num),color="blue"))


ggplot(hyster,aes(as.factor(hysteranthy),mean_size))+geom_boxplot()+stat_summary(color="red")
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")+
  scale_x_discrete(labels=c("seranthous", "hysteranthous"))

ggplot(hyster,aes(as.factor(hysteranthy),mean_num))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flowers/inflorescence")+
  scale_x_discrete(labels=c("seranthous", "hysteranthous"))

ggplot(hyster,aes(as.factor(hysteranthy),display))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("display volume (fl. size*fl. number")+
  scale_x_discrete(labels=c("seranthous", "hysteranthous"))


car::Anova(lm(mean_size~as.factor(hysteranthy),data=hyster),type="III")
car::Anova(lm(mean_num~as.factor(hysteranthy),data=hyster),type="III")
car::Anova(lm(display~as.factor(hysteranthy),data=hyster),type="III")


png("pdsiboxes.png",width = 10,height=10,units = "in",res = 300)
ggplot(d,aes(as.factor(hysteranthy),pdsi))+geom_boxplot()+theme_minimal()+  scale_x_discrete(labels=c("seranthous", "hysteranthous"))
dev.off()

car::Anova(lm(pdsi~as.factor(hysteranthy),data=d),type="III")

?geom_boxplot()
png("pdsimaps.png",width = 10,height=10,units = "in",res = 300)
plot(mean.prunus)
points(lonpoints,latpoints, col=c("royal blue","black"),pch=c(8,5),cex=0.7)
legend(-80,20,c("hysteranthous","seranthous"),pch=c(8,5),col=c("royal blue","black"))
dev.off()

ggplot(d,aes(mean_num,pdsi))+geom_point()+geom_smooth(method="lm",se = TRUE)
ggplot(d,aes(mean_size,pdsi))+geom_point(aes(color=as.factor(hysteranthy)))+geom_smooth(method="lm",se = TRUE)
ggplot(d,aes(display,pdsi))+geom_point()+geom_smooth(method="lm",se = TRUE)
