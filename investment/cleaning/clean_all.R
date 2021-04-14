###prunus data prep and cleaning, see prunus_modeling1.R for other useful code:: April 7 2021 By Dan
# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())



setwd("~/Documents/git/proterant/investment/Input")
d<-read.csv("midwest_round1Dec11.csv") # active datasheet

##subset to useful additions
d.add<-read.csv("species_additions - species_additions.csv")
d.add<-dplyr::filter(d.add,flowers=="Y")
unique(d.add$bbch.f)
##subset to flowering data
d.flo<-dplyr::filter(d,flowers=="Y")

###cleanup mistaken inputs
d.flo<-dplyr::filter(d.flo,bbch.f<80) #clean
d.flo<-dplyr::filter(d.flo,bbch.f>59) #clean
table(d.flo$specificEpithet)
table(d.add$specificEpithet)

d.flo<-rbind(d.flo,d.add) ### this is the data frame with additions

###some 0's were left blank so must be input
d.flo$bbch.v<-ifelse(is.na(d.flo$bbch.v),0,d.flo$bbch.v)

##imput missing values so we dont lose data
table(d.flo$month==0)
table(is.na(d.flo$month))

d.flo$month<-ifelse(d.flo$month==0,4,d.flo$month) ##imputes 12
d.flo$month<-ifelse(is.na(d.flo$month),4,d.flo$month) ## imputes 20

#if day was missing we imput 15

table(d.flo$day==0)
table(is.na(d.flo$day))
d.flo$day<-ifelse(d.flo$day==0,15,d.flo$day) ##imputs 21
d.flo$day<-ifelse(is.na(d.flo$day),15,d.flo$day) #imputes #26

table(d.flo$day>32)## for 3 entries year got input in teh day clolumn...fix it
d.flo$year<-ifelse(d.flo$day>32,d.flo$day,d.flo$year)
d.flo$day<-ifelse(d.flo$day>32,15,d.flo$day)

#clean year and impute 1980 if missing
table(is.na(d.flo$year))
table(d.flo$year==0)

d.flo$year<-ifelse(d.flo$year==2,2002,d.flo$year)
d.flo$year<-ifelse(d.flo$year==5,2005,d.flo$year)
d.flo$year<-ifelse(d.flo$year==6,2006,d.flo$year)
d.flo$year<-ifelse(d.flo$year==198,1980,d.flo$year)
d.flo$year<-ifelse(is.na(d.flo$year),1980,d.flo$year) #imput 27
d.flo$year<-ifelse(d.flo$year==0,1980,d.flo$year) #imput 4
unique(d.flo$year)

###calculate doy from date
d.flo$eventDate2<-paste(d.flo$year,d.flo$month,d.flo$day,sep="-")
d.flo<-filter(d.flo,!is.na(eventDate2))
d.flo$eventDate3<-as.Date(d.flo$eventDate2,format =  "%Y-%m-%d")
d.flo$doy<-yday(d.flo$eventDate3)

###center doy
d.flo$doy.cent<-d.flo$doy-mean(d.flo$doy,na.rm=TRUE)

d.flo<-filter(d.flo,!is.na(doy.cent))
table(d.flo$specificEpithet)

## now rescale BBCH so its even and postive for ordinal regression
unique(d.flo$bbch.v)
d.flo$bbch.v.scale<-NA
d.flo$bbch.v.scale[d.flo$bbch.v==0]<-1
d.flo$bbch.v.scale[d.flo$bbch.v==7]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==9]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==10]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==11]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==14]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==15]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==17]<-5
d.flo$bbch.v.scale[d.flo$bbch.v==19]<-6

##Make FLS only 3 (4) categories to match verbal descriptions
d.flo$bbch.short<-NA
d.flo$bbch.short[d.flo$bbch.v==0]<-1
d.flo$bbch.short[d.flo$bbch.v==7]<-1
d.flo$bbch.short[d.flo$bbch.v==9]<-1
d.flo$bbch.short[d.flo$bbch.v==10]<-2
d.flo$bbch.short[d.flo$bbch.v==11]<-2
d.flo$bbch.short[d.flo$bbch.v==14]<-2
d.flo$bbch.short[d.flo$bbch.v==15]<-2
d.flo$bbch.short[d.flo$bbch.v==17]<-3
d.flo$bbch.short[d.flo$bbch.v==19]<-4


write.csv(d.flo,"input_clean/FLS_clean.csv")

################################
## ADD and clean PDSI data ######
###################################

palmer.b <- brick("..//Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d$lon # make vector of prunus coordinates
latpoints<-d$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
sd.prunus <-calc(palmer.b, fun = sd,na.rm=TRUE) #
min.prunus <-calc(palmer.b, fun = min,na.rm=TRUE) 
ext<-raster::extract(mean.prunus,extract.pts,method="simple")
ext2<-raster::extract(sd.prunus,extract.pts,method="simple")
ext3<-raster::extract(min.prunus,extract.pts,method="simple")
d$pdsi<-ext
d$pdsi.sd<-ext2
d$pdsi.min<-ext3

write.csv(d,"input_clean/pruno_clean_pdsi.csv")

###clean petal length######
petal<-read.csv("flowermeasures - measurements.csv")
petal<-dplyr::select(petal,-X,-X.1,-X.2,-X.3,-X.4)

write.csv(petal,"input_clean/petal_clean.csv")


#### fruit phenology ######
d.fruit<-dplyr::filter(d,fruit=="Y")
table(d.fruit$specificEpithet)

table(d.fruit$month==0)
table(d.fruit$day==0)
table(is.na(d.fruit$month)) ##9
table(is.na(d.fruit$day))
table(is.na(d.fruit$year)) ##8

d.fruit$day<-ifelse(d.fruit$day==0,15,d.fruit$day) ##imputs 5
d.fruit$day<-ifelse(is.na(d.fruit$day),15,d.fruit$day) #imputes #11

d.fruit$year<-ifelse(d.fruit$year==0,2000,d.fruit$year)
d.fruit$year<-ifelse(d.fruit$year==8,2008,d.fruit$year)


###calculate doy from date
d.fruit$eventDate2<-paste(d.fruit$year,d.fruit$month,d.fruit$day,sep="-")
d.fruit<-filter(d.fruit,!is.na(eventDate2))
d.fruit$eventDate3<-as.Date(d.fruit$eventDate2,format =  "%Y-%m-%d")
d.fruit$doy<-yday(d.fruit$eventDate3)

d.fruit<-filter(d.fruit,!is.na(doy))

write.csv(d.fruit,"input_clean/fruit_phen.csv")

###fruit size#######

fruity<-read.csv("fruit_measures - measurements.csv")
fruity<-dplyr::select(fruity,-X,-X.1,-X.2,-X.3,-X.4)

write.csv(fruity,"input_clean/fruitsize_clean.csv")



####freeze tolerance #####
##from https://www.fs.fed.us/rm/boise/AWAE/projects/NFS-regional-climate-change-maps/categories/us-raster-layers.html
str_name<-'..//Data/Tavg_winter_historical.tif' 
winter<-raster(str_name)



extro<-raster::extract(winter,extract.pts,method="simple")
winterT<-(extro-32)*(5/9)
#convert to celcious
d$wintert<-winterT

d %>% group_by(specificEpithet) %>% summarise(meant=mean(wintert,na.rm=TRUE),sdt=sd(wintert,na.rm=TRUE))

write.csv(d,"input_clean/pruno_clean_pdsi_wint.csv")

