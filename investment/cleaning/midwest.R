#### This code 1) reads all prunus data from midwest herbaria
###fixes geo refereneces for as many as possible
### ouputs a sheet of 7838 herbaria sample to score for phenology

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
library(purrr)
library(tidyr)
library(tibble)

mid.herb<-read.csv("SymbOutput_2020-10-26_133412_DwC-A/occurrences.csv")

pruno<-c("alleghaniensis","angustifolia","americana" ,"gracilis","geniculata","hortulana" ,"maritima",
         "mexicana","murrayana","munsoniana","nigra","rivularis","umbellata","subcordata","texana" )

colnames(mid.herb)

mid.herb<-dplyr::select(mid.herb,id,institutionCode,scientificName,specificEpithet,
                        year,month,day,eventDate,country,stateProvince,
                        county,locality,decimalLatitude,decimalLongitude,references)

prunocerasus<-dplyr::filter(mid.herb, specificEpithet %in% pruno)


## sheet 1 with coordinates
pruno.ref<-filter(prunocerasus,!is.na(decimalLatitude))
table(pruno.ref$stateProvince)
table(pruno.ref$specificEpithet)
pruno.ref$lon<-pruno.ref$decimalLongitude
pruno.ref$lat<-pruno.ref$decimalLatitude


#### two without coordinates
pruno.unref<-filter(prunocerasus,is.na(decimalLatitude))

### prep county data
geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" 

colnames(geoCounty)[c(2,7)]<-c("County_name","county")

head(pruno.unref$county)
pruno.unref$county<-tolower(pruno.unref$county) ## make ours love cale



pruno.unref$county<-gsub("county","",pruno.unref$county)
pruno.unref$county<-gsub("()","",pruno.unref$county)
pruno.unref$county<-gsub("co.","",pruno.unref$county,fixed = TRUE)





pruno.unref<-dplyr::left_join(pruno.unref,geoCounty,by=c("county","stateProvince"))
pruno.unref<-dplyr::select(pruno.unref,-County_name,-fips, -state)
intersect(colnames(pruno.unref),colnames(pruno.ref))

checkity<-filter(pruno.unref,!is.na(lon))


pruno.ref$geoclass<-"geo_referenced"
pruno.unref$geoclass<-"county_centroid"

pruneo<-rbind(pruno.ref,pruno.unref)

########## select
#pruneo.geo<-filter(pruneo,!is.na(lat))


 table(pruneo$specificEpithet)


 nested_prun <- pruneo %>%
   group_by(specificEpithet) %>%   # prep for work by Species
   nest() %>%              # --> one row per Species
   ungroup() %>% 
   mutate(n = c(1000,1000,1000,89,1000,457,319,64,682,687,517,246,
                585,188,4))
 
 print(nested_prun)
 sampled_pruny <- nested_prun%>%
   mutate(samp = map2(data, n, sample_n))
 
 sampled_pruny<-sampled_pruny %>% 
   dplyr::select(-data) %>%
   unnest(samp)
 
sampled_pruny<-dplyr::select(sampled_pruny,-n,-locality,-country,-scientificName,-decimalLatitude,-decimalLongitude)

write.csv(sampled_pruny,"midwestherbaria_data_sheet.csv",row.names = FALSE)

######
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
