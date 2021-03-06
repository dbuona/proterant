### Started by Dan Feb 9 2021
### Part 1:### modeling FLS for FNA with measurement models in brms
#FLS ~ pdsi + flower traits + fruit traits
###Part 2: Focusing on prunocerasus with stan model
#based on 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library("housingData")
library(stringr)
library("ncdf4")
library(raster)
library(ggrepel)
library(patchwork)
library(brms)
library("grid")
library(tidybayes)



#-------------------------#
#------Part 1--------------#
#--------------------------#
setwd("~/Documents/git/proterant/investment")
#load("FNA.Rda")
FNA<-read.csv("Data/cherry_data.csv") ## measurement and FLS data from FNA
##clean specicies name

FNA$species<-ifelse(FNA$species=="hortunlana","hortulana",FNA$species)
FNA$species<-ifelse(FNA$species=="gladulosa","glandulosa",FNA$species)
FNA$species<-ifelse(FNA$species=="speciosa","speciosa",FNA$species)
FNA$species<-ifelse(FNA$species=="fasiculata","fasciculata",FNA$species)

#### get herbaria specimen coordinates
mid.herb<-read.csv("SymbOutput_2020-10-26_133412_DwC-A/occurrences.csv")
mid.herb<-filter(mid.herb,specificEpithet %in% unique(FNA$species)) ## filter to Species we have data 4
unique(mid.herb$specificEpithet)
## select rows that already have coordinate
pruno.ref<-filter(mid.herb,!is.na(decimalLatitude))
table(pruno.ref$stateProvince)
table(pruno.ref$specificEpithet)
pruno.ref$lon<-pruno.ref$decimalLongitude
pruno.ref$lat<-pruno.ref$decimalLatitude

##now add county level coordiantes
#### two without coordinates
pruno.unref<-filter(mid.herb,is.na(decimalLatitude)) ## filter entries with no lat/lon

### prep county data
geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" 

colnames(geoCounty)[c(2,7)]<-c("County_name","county") ## MATCH NAMES

pruno.unref<-filter(pruno.unref,id!=17453277) #remvoe problematic canadian entry
pruno.unref$county<-tolower(pruno.unref$county) ## make ours lower case



pruno.unref$county<-gsub("county","",pruno.unref$county) ### get rid of extenious coounty info
pruno.unref$county<-gsub("()","",pruno.unref$county) # ""
pruno.unref$county<-gsub("co.","",pruno.unref$county,fixed = TRUE) # ""




pruno.unref<-dplyr::left_join(pruno.unref,geoCounty,by=c("county","stateProvince"))
pruno.unref<-dplyr::select(pruno.unref,-County_name,-fips, -state)
intersect(colnames(pruno.unref),colnames(pruno.ref))

pruno.ref$geoclass<-"geo_referenced"
pruno.unref$geoclass<-"county_centroid"
pruno.unref<-filter(pruno.unref,!is.na(lon)) # select only entries with lat lon
pruneo<-rbind(pruno.ref,pruno.unref)

##now add pdsi data
palmer.b <- brick("Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-pruneo$lon # make vector of prunus coordinates
latpoints<-pruneo$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
sd.prunus <-calc(palmer.b, fun = sd,na.rm=TRUE) #
min.prunus <-calc(palmer.b, fun = min,na.rm=TRUE) 
ext<-raster::extract(mean.prunus,extract.pts,method="simple")
ext2<-raster::extract(sd.prunus,extract.pts,method="simple")
ext3<-raster::extract(min.prunus,extract.pts,method="simple")
pruneo$pdsi<-ext
pruneo$pdsi.sd<-ext2
pruneo$pdsi.min<-ext3
###now do winter T
str_name<-'Data/Tavg_winter_historical.tif' 
winter<-raster(str_name)



extro<-raster::extract(winter,extract.pts,method="simple")
pruneo$winterT<-(extro-32)*(5/9)

#pdsi summary data
minT<- pruneo %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanT=mean(winterT,na.rm=TRUE),sdT=sd(winterT,na.rm=TRUE))
colnames(minT)[1]<-"species"
setdiff(FNA$species,pdsi$species) ## lose speciosa

pdsi<- pruneo %>% dplyr::group_by(specificEpithet) %>% dplyr::summarise(meanpdsi=mean(pdsi,na.rm=TRUE),sdmean=sd(pdsi,na.rm=TRUE),minpdsi=mean(pdsi.min,na.rm=TRUE),sdmin=sd(pdsi.min,na.rm=TRUE))
colnames(pdsi)[1]<-"species"
setdiff(FNA$species,pdsi$species) ## lose speciosa


## combine data
FNA<-left_join(FNA,pdsi)
FNA<-left_join(FNA,minT)
##quick plots
cor(FNA$petal_low,FNA$petal_high) # .88
cor(FNA$inflor_low,FNA$inflor_high) #.86 can probably use jsut one or the other in the
cor(FNA$inflor_low,FNA$petal_low) #-.35 right direction but not as correlated as trade off would suggest
cor(FNA$inflor_high,FNA$petal_high) #-.31
cor(FNA$fruit_low,FNA$fruit_high) #.88
##Make FLS numeric for oridnal modeling
unique(FNA$FLS)
FNA$FLSnum<-NA
FNA$FLSnum[FNA$FLS=="before"]<-1
FNA$FLSnum[FNA$FLS=="before/with"]<-2
FNA$FLSnum[FNA$FLS=="with"]<-3
FNA$FLSnum[FNA$FLS=="after"]<-4
FNA$FLSnum<-as.integer(FNA$FLSnum)
##quick plots
class(FNA$FLSnum)
plot(FNA$FLSnum,FNA$inflor_low)

##basic model
# FLSnum~ inflorlow*petal_low+me(pdsi,sd)+fruit_low +phylo someday

#rough sd and mean calculation for other predictons
#FNA$petalmean<-(FNA$petal_low+FNA$petal_high)/2
#FNA$petalsd<-(FNA$petal_high-FNA$petal_low)/4 #"rule of ranges"


#FNA$inflormean<-(FNA$inflor_low+FNA$inflor_high)/2
#FNA$inflorsd<-(FNA$inflor_high-FNA$inflor_low)/4 #"rule of ranges"

#FNA$fruitmean<-(FNA$fruit_low+FNA$fruit_high)/2
#FNA$fruitsd<-(FNA$fruit_high-FNA$fruit_low)/4 #"rule of ranges"

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
#FNA$fruit.z.mean<-range01(FNA$fruitmean)
#FNA$fruit.z.sd<-range01(FNA$fruitsd)

#FNA$inflor.z.mean<-range01(FNA$inflormean)
#FNA$inflor.z.sd<-range01(FNA$inflorsd)

#FNA$petal.z.mean<-range01(FNA$petalmean)
#FNA$petal.z.sd<-range01(FNA$petalsd)

#zscore all predictor
FNA$meanpdsi.z<-zscore(FNA$meanpdsi)
FNA$minpdsi.z<-zscore(FNA$minpdsi)
FNA$petal.z<-zscore(FNA$petal_high)
FNA$inflor.z<-zscore(FNA$inflor_high)
FNA$fruit.z<-zscore(FNA$fruit_high)
FNA$cold.z<-zscore(FNA$meanT)


#get_prior(FLSnum~me(petalmean,petalsd)+me(meanpdsi,sdmean),data=FNA)


FNAordz<-brm(FLSnum~petal.z+inflor.z+fruit.z+cold.z+meanpdsi.z,
             data=FNA,
             family=cumulative("logit"),warmup=3000,iter=4000)

fixef(FNAordz)
#FNAord.traits<-brm(FLSnum~me(petalmean,petalsd)+me(inflormean,inflorsd)+
 #                    me(fruitmean,fruitsd),
  #                 data=FNA,
   #                family=cumulative("logit"),warmup=3000,iter=4000,control = list(adapt_delta=.95))

launch_shinystan(FNAordz.error)


#get_prior(FLSnum~petal.z.mean+inflor.z.mean+
 # me(meanpdsi.z,sdpdsi.z)+fruit.z.mean,
#data=FNA)

#bprior <- c(prior(normal(0,1), class = b),
 #           prior(normal(0,1), class = meanme ),
  #          prior(uniform(0,1), class = sdme), 
   #               prior(normal(0,1), class = Intercept))

#FNAordz.error1<-brm(FLSnum~petal.z.mean+inflor.z.mean+
 #                    meanpdsi.z+me(fruit.z.mean,fruit.z.sd),
  #                 data=FNA ,save_mevars = TRUE,
   #                family=cumulative("logit"),warmup=3000,iter=4000 )

#launch_shinystan(FNAordz.error1)
p1<-plot(conditional_effects(FNAordz, "meanpdsi.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p5<-plot(conditional_effects(FNAordz, "cold.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p2<-plot(conditional_effects(FNAordz, "petal.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p3<-plot(conditional_effects(FNAordz, "inflor.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p4<-plot(conditional_effects(FNAordz, "fruit.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))

p1<-p1[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p2<-p2[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p3<-p3[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p4<-p4[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p5<-p5[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))


jpeg("Plots/FNA_mean_ordinal.jpeg")
ggpubr::ggarrange(p1,p5,p2,p3,p4,common.legend = TRUE)
dev.off()
#### run this just on prunocerasus
pruno<-c("alleghaniensis","angustifolia","americana" ,"gracilis","geniculata","hortulana" ,"maritima",
         "mexicana","murrayana","munsoniana","nigra","rivularis","umbellata","subcordata","texana" )


FNA.pruno<-dplyr::filter(FNA, species %in% pruno)

FNAordz.pruno<-brm(FLSnum~petal.z+fruit.z+cold.z+meanpdsi.z,
             data=FNA.pruno,
             family=cumulative("logit"),warmup=3000,iter=4000)

fixef(FNAordz.pruno)


p11<-plot(conditional_effects(FNAordz.pruno, "meanpdsi.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p55<-plot(conditional_effects(FNAordz.pruno, "cold.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p22<-plot(conditional_effects(FNAordz.pruno, "petal.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p44<-plot(conditional_effects(FNAordz.pruno, "fruit.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))

p11<-p11[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p22<-p22[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p44<-p44[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p55<-p55[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))

jpeg("Plots/FNA_mean_prunocerasus.jpeg")
ggpubr::ggarrange(p11,p55,p22,p44,common.legend = TRUE)
dev.off()
fixef(FNAordz.pruno,probs = c(.25,.75))

cor(FNA$petal_low,FNA$meanpdsi,use = "complete.obs")
cor(FNA$minpdsi,FNA$meanpdsi,use = "complete.obs")

cor(FNA.pruno$minpdsi,FNA.pruno$meanpdsi,use = "complete.obs")
cor(FNA.pruno$meanT,FNA.pruno$meanpdsi,use = "complete.obs")

## make sure we can use the Zanne tree
colnames(FNA)
FNA$name<-paste("Prunus",FNA$species, sep="_")

library(ape)
library(phytools)
library(geiger)
library(caper)
library(picante)
library(tidyverse)
library(boot)
library(phylolm)



treee<-read.tree("..//input/Vascular_Plants_rooted.dated.tre")
names.intree<-treee$tip.label
namelist<-FNA$name

which(names.intree%in%namelist)
length(which(namelist%in%names.intree))
length(which(!namelist%in%names.intree))
save.image("FNA.Rda")


