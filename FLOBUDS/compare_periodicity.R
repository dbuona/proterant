###compare_periodicity begun by Dan B. on Feb 7 2019
#The goal of this file is to compare the effect of covarying photoperiodicity and thermoperiodcity by comparing dany flynn to Dan b's budburst experients

##Ask Lizzie about the chilling.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(brms)

setwd("~/Documents/git/proterant/FLOBUDS")

db.dat<-read.csv("first.event.dat.csv")
df.dat<-read.csv("input/Budburst By Day.csv")

df.dat.HF<-filter(df.dat,site=="HF") ##filter to only HF species
n <- 4
df.dat.HF$GEN.SPA<-paste(substr(df.dat.HF$sp, 1, n-1), ".", substr(df.dat.HF$sp, n, nchar(df.dat.HF$sp)), sep = "") ## make a column for GEN.SA


df.dat.HF<-filter(df.dat.HF,chill=="chill0")### select only DF's no chill treatment assume there was field chillint


db.dat$photoperiod<-ifelse(db.dat$Light=="L",12,8)
db.dat$temp_day<-ifelse(db.dat$Force=="W",24,18)
db.dat$temp_night<-ifelse(db.dat$Force=="W",18,12)
db.dat$chilldays<-ifelse(db.dat$Chill==0,28,56)

unique(db.dat$GEN.SPA)
unique(df.dat.HF$GEN.SPA)
matching_sps<-intersect(unique(df.dat.HF$GEN.SPA),unique(db.dat$GEN.SPA)) ## could also add viburnum and vaccinium at the genus level


df.dat.match<-df.dat.HF %>%filter(GEN.SPA %in% c(matching_sps))  ### make both dataset noly consiston matching sps
db.dat.match<-db.dat %>%filter(GEN.SPA %in% c(matching_sps))





 ### select only DAn B's low chilling (28days)
### Dan Flynn Temps  20/10 (15) and 15/5  (10)
##Dan B's Temps: 24/18 (21) and 18/12 (15)

##Dan B and Dan F photoperiods 12hours and 8 hours. Dan F covaried theromo anf photoperiodicity, Dan B did not.
###Chilling. Dan Flyyn's "Chill 2" was 4 degrees according to paper 33 days at 4c
##Dan Flynn's chilling periods were 0 days 33 days at 4, adn 33 days at 2.

#ubset DF's data to be similar to DB's


db.dat.match$FORCE<-ifelse(db.dat.match$temp_day==24,1,0) ### new column where forcing temperature are binary low or high
df.dat.match$FORCE<-ifelse(df.dat.match$warm==20,1,0)


db.dat.match$study<-"DB"
df.dat.match$study<-"DF"

##subset to columns of use
db.dat.match<-dplyr::select(db.dat.match,GEN.SPA,flo_day.60.,Lbb_day.9.,leaf_day.15.,treatcode,temp_day,photoperiod,Chill,study,FORCE)
df.dat.match<-dplyr::select(df.dat.match,GEN.SPA,fday,bday,lday,treatcode,warm,photo,chill,study,FORCE)




###make columns the same
colnames(db.dat.match)[colnames(db.dat.match)=="flo_day.60." ] <- "fday"
colnames(db.dat.match)[colnames(db.dat.match)=="leaf_day.15." ] <- "lday"
colnames(db.dat.match)[colnames(db.dat.match)=="Lbb_day.9." ] <- "bday"
colnames(db.dat.match)[colnames(db.dat.match)=="temp_day" ] <- "warm"
colnames(db.dat.match)[colnames(db.dat.match)=="photoperiod" ] <- "photo"
colnames(db.dat.match)[colnames(db.dat.match)=="Chill" ] <- "chill"

colnames(db.dat.match)
colnames(df.dat.match)

unique(db.dat.match$chill)
db.dat.match.short<-filter(db.dat.match,chill==0)
db.dat.match.long<-filter(db.dat.match,chill==1)

both.dat1<-rbind(df.dat.match,db.dat.match.short) ###This is a datasheet
both.dat2<-rbind(df.dat.match,db.dat.match.long)

both.dat1$PHOTO<-ifelse(both.dat1$photo==8,0,1)
both.dat2$PHOTO<-ifelse(both.dat2$photo==8,0,1)

both.dat1$warm.cent<-both.dat1$warm-mean(both.dat1$warm)
both.dat1$photo.cent<-both.dat1$photo-mean(both.dat1$photo)
both.dat2$warm.cent<-both.dat2$warm-mean(both.dat2$warm)
both.dat2$photo.cent<-both.dat2$photo-mean(both.dat2$photo)


both.mod.bud1<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat1)
summary(both.mod.bud1)

both.mod.bud2<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat2)

summary(both.mod.bud2)

both.mod.leaf1<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat1)
summary(both.mod.leaf1)

both.mod.leaf2<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat2)
summary(both.mod.leaf2)

both.mod.bud1.int<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat1)
summary(both.mod.bud1.int)

both.mod.leaf1.int<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat1)
summary(both.mod.leaf1.int)

both.mod.bud2.int<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat2)
summary(both.mod.bud2.int)

both.mod.leaf2.int<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat2)
summary(both.mod.leaf2.int)

####try it with continuous preidctors
both.mod.bud1.cont<-brm(bday~study+warm.cent+photo.cent+study:warm.cent+study:photo.cent+(study+warm.cent+photo.cent+study:warm.cent+study:photo.cent|GEN.SPA),data=both.dat1)
summary(both.mod.bud1.cont)






####but above effects could jsut be cause Dan B used consitantly higher forcing. Luckily my low force average is teh same as DF's high average


both.dat.sameforce2<-both.dat2 %>% filter(warm %in% c(20,18))# This make a data sheet where the average forcing in both experiments was 15 C (DB low, DF high)
both.dat.sameforce1<-both.dat1 %>% filter(warm %in% c(20,18))

photo.only.bb1<-brm(bday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2000)
summary(photo.only.bb1)
coef(photo.only.bb1)

photo.only.leaf1<-brm(lday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2000)
summary(photo.only.leaf1)

photo.only.bb2<-brm(bday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce2,iter=3000,warmup=2000)
summary(photo.only.bb2)
coef(photo.only.bb2)

photo.only.leaf2<-brm(lday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce2,iter=3000,warmup=2000)
summary(photo.only.leaf2)




