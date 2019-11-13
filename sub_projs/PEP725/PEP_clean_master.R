###PEP cleaning file
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/proterant/sub_projs/PEP725")

library(dplyr)
library(tidyr)
###seperate the one column into useful ones

###function for formating data
colfix<-function(x){names(x)[names(x)=="s_id.lon.lat.alt.plant_id.cult_id.bbch.year.day"] <- "s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day" #change names
x<-separate(x,"s_id;lon;lat;alt;plant_id;cult_id;bbch;year;day",into= c("s_id","lon","lat","alt","plant_id","cult_id","bbch","year","day"), sep=";" ) # seperate
}

makedaters<-function(x,y){
  z<-rbind(x,y)
  z<-unite(z,locale,lat,lon,sep=",",remove=FALSE)
  z$bbch<-ifelse(z$bbch=="11","leaf_day","flo_day")
  z<-spread(z,bbch,day)
  z<-transform(z,flo_day=as.numeric(flo_day))
  z<-transform(z,leaf_day=as.numeric(leaf_day))
  z<-transform(z,FLS=leaf_day-flo_day)
  z<-filter(z,FLS<50)
  z<-filter(z,FLS>(-70))
}

siteagg50 <-function(x){y<- aggregate(x[("year")], x[c("s_id")],FUN=length) #count how many year of observartion there are at each site
y<-subset(y, year>50)}

siteagg20 <-function(x){y<- aggregate(x[("year")], x[c("s_id")],FUN=length) #count how many year of observartion there are at each site
y<-subset(y, year>20)}

datify<-function(x,y){x[which(x$s_id %in% y$s_id),]}
########SPECIES######################################################################

###Fraxinus
frax.leaf<-read.csv("pep_raw/buo_120_000_011.csv", header=TRUE)
frax.flo<-read.csv("pep_raw/buo_120_000_060.csv",header=TRUE)

frax.leaf<-colfix(frax.leaf) 
frax.flo<-colfix(frax.flo)
frax.dat<-makedaters(frax.leaf,frax.flo)

frax50 <- siteagg50(frax.dat)
frax20 <- siteagg20(frax.dat)

frax.dat.20<- datify(frax.dat,frax20)
frax.dat.50<-datify(frax.dat,frax50)


########Aesculus
aesc.leaf<-read.csv("pep_raw/buo_101_000_011.csv", header=TRUE)
aesc.flo<-read.csv("pep_raw/buo_101_000_060.csv",header=TRUE)

aesc.leaf<-colfix(aesc.leaf)
aesc.flo<-colfix(aesc.flo)
aesc.dat<-makedaters(aesc.leaf,aesc.flo)

aesc50 <- siteagg50(aesc.dat)
aesc20 <- siteagg20(aesc.dat)

aesc.dat.20<- datify(aesc.dat,aesc20)
aesc.dat.50<- datify(aesc.dat,aesc50)

##Alnus glutinosa
alnu.leaf<-read.csv("pep_raw/buo_102_040_011.csv", header=TRUE)
alnu.flo<-read.csv("pep_raw/buo_102_040_060.csv",header=TRUE)

alnu.leaf<-colfix(alnu.leaf)
alnu.flo<-colfix(alnu.flo)
alnu.dat<-makedaters(alnu.leaf,alnu.flo)

alnu50 <- siteagg50(alnu.dat)
alnu20 <- siteagg20(alnu.dat)

alnu.dat.20<- datify(alnu.dat,alnu20)
alnu.dat.50<- datify(alnu.dat,alnu50)


##FAGUS
fagu.leaf<-read.csv("pep_raw/buo_108_010_011.csv", header=TRUE)
fagu.flo<-read.csv("pep_raw/buo_108_010_060.csv",header=TRUE)

fagu.leaf<-colfix(fagu.leaf)
fagu.flo<-colfix(fagu.flo)
fagu.dat<-makedaters(fagu.flo,fagu.leaf)

fagu50 <- siteagg50(fagu.dat)## not enough for usage
fagu20 <- siteagg20(fagu.dat)

fagu.dat.20<- datify(fagu.dat,fagu20)


#Tilia
tili.leaf<-read.csv("pep_raw/buo_129_070_011.csv", header=TRUE)
tili.flo<-read.csv("pep_raw/buo_129_070_060.csv",header=TRUE)

tili.leaf<-colfix(tili.leaf)
tili.flo<-colfix(tili.flo)
tili.dat<-makedaters(tili.flo,tili.leaf)

tili50 <- siteagg50(tili.dat)## not enough for usage
tili20 <- siteagg20(tili.dat)

tili.dat.20<- datify(tili.dat,tili20)

###Corylus
cory.leaf<-read.csv("pep_raw/buo_107_000_011.csv", header=TRUE)
cory.flo<-read.csv("pep_raw/buo_107_000_060.csv",header=TRUE)

cory.leaf<-colfix(cory.leaf)
cory.flo<-colfix(cory.flo)
cory.dat<-makedaters(cory.leaf,cory.flo)

cory50 <- siteagg50(cory.dat)## no enough for usage
cory20 <- siteagg20(cory.dat) ## nodata for use


###Betula pendula
betu.leaf<-read.csv("pep_raw/buo_106_020_011.csv", header=TRUE)
betu.flo<-read.csv("pep_raw/buo_106_020_060.csv",header=TRUE)

betu.leaf<-colfix(betu.leaf)
betu.flo<-colfix(betu.flo)
betu.dat<-makedaters(betu.leaf,betu.flo)

betu50 <- siteagg50(betu.dat)## no enough for usage
betu20 <- siteagg20(betu.dat) ## nodata for use

betu.dat.20<- datify(betu.dat,betu20)
betu.dat.50<- datify(betu.dat,betu50)

##Acer platanus
acer.leaf<-read.csv("pep_raw/buo_115_030_011.csv", header=TRUE)
acer.flo<-read.csv("pep_raw/buo_115_030_060.csv",header=TRUE)

acer.leaf<-colfix(acer.leaf)
acer.flo<-colfix(acer.flo)
acer.dat<-makedaters(acer.leaf,acer.flo)

acer50 <- siteagg50(acer.dat)## no enough for usage
acer20 <- siteagg20(acer.dat) ##none

## write data sheets#######

write.csv(frax.dat.20,"input/frax20year.csv")
write.csv(frax.dat.50,"input/frax50year.csv")

write.csv(aesc.dat.20,"input/aesc20year.csv")
write.csv(aesc.dat.50,"input/aesc50year.csv")

write.csv(alnu.dat.20,"input/alnu20year.csv")
write.csv(alnu.dat.50,"input/alnu50year.csv")

write.csv(fagu.dat.20,"input/fagu20year.csv")

write.csv(tili.dat.20,"input/tili20year.csv")

write.csv(betu.dat.20,"input/betu20year.csv")
write.csv(betu.dat.50,"input/betu50year.csv")
