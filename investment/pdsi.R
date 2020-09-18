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


###hysteranthy
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

table(hyster$FLS)
table(hyster$hyst)
table(hyster$hyst2)
prunus.data.2<-dplyr::left_join(prunus.data,hyster,by="specificEpithet")


d<-dplyr::filter(prunus.data.2,!is.na(FLS)) ### 7668 rows have county coordiates






## hypothesis1: hysteranthy relates to flower size

ggplot(hyster,aes(inflor_type))+geom_bar(position="dodge",aes(fill=hyst))


a<-ggplot(hyster,aes(FLS,mean_size))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")

aa<-ggplot(hyster,aes(hyst,mean_size))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flower size (mm)")


b<-ggplot(hyster,aes(FLS,mean_num))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flowers/inflorescence")
bb<-ggplot(hyster,aes(hyst,mean_num))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean flowers/inflorescence")


c<-ggplot(hyster,aes(FLS,display))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("display volume (fl. size*fl. number")
cc<-ggplot(hyster,aes(hyst,display))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("display volume (fl. size*fl. number)")


e<-ggplot(hyster,aes(FLS,mean_fruit))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean fruit siam (mm)")
ee<-ggplot(hyster,aes(hyst,mean_fruit))+geom_boxplot()+stat_summary(color="red")+
  ggthemes::theme_base()+xlab("FLS")+ylab("mean fruit size (mm)")



ggpubr::ggarrange(a,b,c,e)
ggpubr::ggarrange(aa,bb,cc,ee)




mod1<-brms::brm(mean_fruit~hyst,data=hyster)


0car::Anova(lm(mean_size~FLS,data=hyster),type="III")


car::Anova(lm(mean_num~FLS,data=hyster),type="III")
car::Anova(lm(mean_num~hyst,data=hyster),type="III")


car::Anova(lm(display~FLS,data=hyster),type="III")
car::Anova(lm(display~hyst,data=hyster),type="III")

car::Anova(lm(mean_fruit~FLS,data=hyster),type="III")
car::Anova(lm(mean_fruit~hyst,data=hyster),type="III")


png("pdsiboxes.png",width = 10,height=10,units = "in",res = 300)
ggplot(d,aes(hyst,pdsi))+geom_boxplot()+theme_minimal()

car::Anova(lm(pdsi~hyst,data=d),type="III")

?geom_boxplot()
png("pdsimaps.png",width = 10,height=10,units = "in",res = 300)
plot(mean.prunus)
points(lonpoints,latpoints, col=c("royal blue","black"),pch=c(8,5),cex=0.7)
legend(-80,20,c("hysteranthous","seranthous"),pch=c(8,5),col=c("royal blue","black"))
dev.off()

ggplot(d,aes(mean_num,pdsi))+geom_point(aes(color=FLS))+geom_smooth(method="lm",se = TRUE)
ggplot(d,aes(mean_size,pdsi))+geom_point(aes(color=hyst))+geom_smooth(method="lm",se = TRUE)
ggplot(d,aes(display,pdsi))+geom_point(aes(color=hyst))+geom_smooth(method="lm",se = TRUE)
