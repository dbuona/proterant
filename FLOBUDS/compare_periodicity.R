###compare_periodicity begun by Dan B. on Feb 7 2019
#The goal of this file is to compare the effect of covarying photoperiodicity and thermoperiodcity by comparing dany flynn to Dan b's budburst experients

##Ask Lizzie about the chilling.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(brms)

setwd("~/Documents/git/proterant/FLOBUDS")

db.dat<-read.csv("input/xfirst.event.dat.csv")
df.dat<-read.csv("input/hf314-01-budburst.csv")

df.dat.HF<-filter(df.dat,site=="HF")
n <- 4
df.dat.HF$GEN.SPA<-paste(substr(df.dat.HF$sp, 1, n-1), ".", substr(df.dat.HF$sp, n, nchar(df.dat.HF$sp)), sep = "")
df.dat.HF<-filter(df.dat.HF,chill=="chill2")


unique(db.dat$GEN.SPA)
unique(df.dat.HF$GEN.SPA)
matching_sps<-intersect(unique(df.dat.HF$GEN.SPA),unique(db.dat$GEN.SPA)) ##



db.dat$photoperiod<-ifelse(db.dat$Light=="L",12,8)
db.dat$temp_day<-ifelse(db.dat$Force=="W",24,18)
db.dat$temp_night<-ifelse(db.dat$Force=="W",18,12)
db.dat$chilldays<-ifelse(db.dat$Chill==0,28,56)


### Dan Flynn Temps  20/10 (15) and 15/5  (10)
##Dan B's Temps: 24/18 (21) and 18/12 (15)

##Dan B and Dan F photoperiods 12hours and 8 hours. Dan F covaried theromo anf photoperiodicity, Dan B did not.

###Chilling. Dan Flyyn's "Chill 2" was 4 degrees according to paper 33 days at 4c
##Dan Flynn's chilling periods were 0, 30d and 45d, Dan B's were 28d, and 56d. Acutally I think chill 2 was 4 degrees at

#ubset DF's data to be similar to DB's
df.dat<-filter(df.dat,chill=="chill2")

db.dat$chill[db.dat$Chill==0]<-"chill1"
db.dat$chill[db.dat$Chill==1]<-"chill2"


### try with HF only
df.dat.HF<-filter(df.dat,site=="HF")
unique(db.dat$GEN.SPA)
unique(df.dat.HF$sp)

n <- 4
df.dat.HF$GEN.SPA<-paste(substr(df.dat.HF$sp, 1, n-1), ".", substr(df.dat.HF$sp, n, nchar(df.dat.HF$sp)), sep = "")


###use only overlapping species
 ### overlapping speices

db.dat.match<-filter(db.dat,GEN.SPA %in% matching_sps)
df.dat.hf.match<-filter(df.dat.HF,GEN.SPA %in% matching_sps)

db.dat.match$study<-"DB"
df.dat.hf.match$study<-"DF"

db.dat.match<-dplyr::select(db.dat.match,GEN.SPA,flo_day,Lbb_day,leaf_day,treatcode,photoperiod,temp_day,temp_night, chilldays,chill,study)

df.dat.hf.match<-dplyr::select(df.dat.hf.match,GEN.SPA,fday,bday,lday,treatcode,warm,photo,chill,study)

df.dat.hf.match$chilldays<-33
#ifelse(df.dat.hf.match$chill=="chill1",30,45) ignore this
df.dat.hf.match$temp_night<-ifelse(df.dat.hf.match$warm==20,10,5)


###make columns the same


colnames(db.dat.match)[colnames(db.dat.match)=="flo_day" ] <- "fday"
colnames(db.dat.match)[colnames(db.dat.match)=="leaf_day" ] <- "lday"
colnames(db.dat.match)[colnames(db.dat.match)=="Lbb_day" ] <- "bday"
colnames(db.dat.match)[colnames(db.dat.match)=="temp_day" ] <- "warm"
colnames(db.dat.match)[colnames(db.dat.match)=="photoperiod" ] <- "photo"

### make db's dat a relative matching chill period
db.dat.match<-filter(db.dat.match,chilldays==28)

both.dat<-rbind(df.dat.hf.match,db.dat.match)

both.mod<-brm(bday~warm+photo+warm:study+photo:study+warm:ph(1|GEN.SPA),data=both.dat)
coef(both.mod)
summary(both.mod)
