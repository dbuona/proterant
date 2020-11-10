rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()



setwd("~/Documents/git/proterant/investment")
library("housingData")
library(stringr)
library("ncdf4")
library(raster)
library(ggplot2)
library(dplyr)
mid.herb<-read.csv("SymbOutput_2020-10-26_133412_DwC-A/occurrences.csv")

pruno<-c("alleghaniensis","angustifolia","americana" ,"gracilis","geniculata","hortulana" ,"maritima",
         "mexicana","murrayana","munsoniana","nigra","rivularis","umbellata","subcordata","texana" )

colnames(mid.herb)
prunocerasus<-dplyr::filter(mid.herb, specificEpithet %in% pruno)

pruno.ref<-filter(prunocerasus,!is.na(decimalLatitude))
pruno.unref<-filter(prunocerasus,is.na(decimalLatitude))

pruno.ref$lon<-pruno.ref$decimalLongitude
pruno.ref$lat<-pruno.ref$decimalLatitude


geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties

colnames(geoCounty)[6]<-"stateProvince" 

head(geoCounty)
head(pruno.unref$county)
pruno.unref$county<-paste(pruno.unref$county,"County",sep=" ")


pruno.unref<-dplyr::left_join(pruno.unref,geoCounty,by=c("county","stateProvince"))
pruno.unref<-dplyr::select(pruno.unref,-rMapCounty,-fips, -state)
colnames(pruno.unref)

pruneo<-rbind(pruno.ref,pruno.unref)


palmer.b <- brick("lbda-v2_kddm_pmdi_2017.nc")  ## read in balmer drought index data 

lonpoints<-pruneo$lon # make vector of prunus coordinates
latpoints<-pruneo$lat #


extract.pts <- cbind(lonpoints,latpoints) #make coordinates to extract

mean.prunus <- calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
plot(mean.prunus)
ext <- raster::extract(mean.prunus,extract.pts,method="simple")
 
pruneo$pdsi<-ext
table(pruneo$institutionCode)

hyster<-read.csv("Data/cherry_data.csv")

colnames(hyster)[1]<-"specificEpithet"
unique(hyster$specificEpithet)

hyster$specificEpithet<-ifelse(hyster$specificEpithet=="hortunlana","hortulana",hyster$specificEpithet)

#pruneo$specificEpithet<-ifelse(pruneo$specificEpithet=="munsoniana","rivularis",pruneo$specificEpithet)
#pruneo$specificEpithet<-ifelse(pruneo$specificEpithet=="alleghaniensis","umbellata",pruneo$specificEpithet)

pruneo2<-left_join(pruneo,hyster)

a<-ggpubr::ggboxplot('specificEpithet','pdsi',data=pruneo2,color='FLS')
b<-ggpubr::ggboxplot('FLS','pdsi',data=pruneo2,color='FLS')
mod<-lm(pdsi~FLS,data=pruneo2)
car::Anova(mod,type="III")

HSD.test(mod,"FLS",group=TRUE,console=TRUE)
### how many speciments
specs<-geo %>%group_by(specificEpithet)%>% count() %>% arrange(desc(n))
sum(specs$n)

before<-c("texana","geniculata","umbellata","gracilis", "angustifolia","nigra","maritima","alleghaniensis")
bw<-c("rivularis","mexicana","munsoniana")
with<-c("subcordata","hortulana","americana","murrayana")

pruneo2$FLS2<-NA

pruneo2$FLS2[which(pruneo2$specificEpithet %in% before)]<-"before"
pruneo2$FLS2[which(pruneo2$specificEpithet %in% bw)]<-"before/with"
pruneo2$FLS2[which(pruneo2$specificEpithet %in% with)]<-"with"
pruneo2$FLS2<-as.factor(pruneo2$FLS2)

aa<-ggpubr::ggboxplot('specificEpithet','pdsi',data=pruneo2,color='FLS2')
bb<-ggpubr::ggboxplot('FLS2','pdsi',data=pruneo2,color='FLS2')

mod<-lm(pdsi~FLS2,data=pruneo2)
car::Anova(mod,type="III")
aHSD.test(mod,"FLS2",group=TRUE,console=TRUE)

ggpubr::ggarrange(a,aa,nrow=2)
