###compare_periodicity begun by Dan B. on Feb 7 2019
#The goal of this file is to compare the effect of covarying photoperiodicity and thermoperiodcity by comparing dany flynn to Dan b's budburst experients

##Ask Lizzie about the chilling.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(brms)
library(chillR)

setwd("~/Documents/git/proterant/FLOBUDS")

db.dat<-read.csv("Flo_buds_for_thermoperiodicity.csv")
df.dat<-read.csv("input/Budburst By Day.csv")
##
##what chilling to use  in my experiments
#dayseq = seq(as.numeric(format(as.POSIXlt("2017-10-26", "%Y-%m-%d"), "%j")),
 #            as.numeric(format(as.POSIXlt("2017-11-21", "%Y-%m-%d"), "%j")))

#chillDB <- data.frame(
#  Year = as.numeric(rep("2017", 27*24)),
 # JDay = as.numeric(rep(dayseq, each = 24)),
  #Hour = rep(0:23, 27),
  #Temp = 4
#)



#dayseq2 = seq(as.numeric(format(as.POSIXlt("2017-10-26", "%Y-%m-%d"), "%j")),
 #            as.numeric(format(as.POSIXlt("2017-12-19", "%Y-%m-%d"), "%j")))

#chill2 <- data.frame(
 #Year = as.numeric(rep("2017", 55*24)),
#  JDay = as.numeric(rep(dayseq2, each = 24)),
 # Hour = rep(0:23, 55),
  #Temp = 4
#)


# bug: Jday vs JDay, returns error for "chillout not found"


#chill1 <- chilling(chillDB, 300, 326)
#chill.2 <- chilling(chill2, 300, 354)


#DFfield 819 chillhours
#DB1: 624 chill hours
#DB2: 1,296 chill hour

df.dat.HF<-filter(df.dat,site=="HF") ##filter to only HF species
n <- 4
df.dat.HF$GEN.SPA<-paste(substr(df.dat.HF$sp, 1, n-1), ".", substr(df.dat.HF$sp, n, nchar(df.dat.HF$sp)), sep = "") ## make a column for GEN.SA


df.dat.HF<-filter(df.dat.HF,chill=="chill0")### select only DF's no chill treatment assume there was field chillint

db.dat$photoperiod<-ifelse(db.dat$Light=="L",12,8)
db.dat$temp_day<-ifelse(db.dat$Force=="W",24,18)
db.dat$temp_night<-ifelse(db.dat$Force=="W",18,12)
db.dat$chilldays<-ifelse(db.dat$Chill==0,28,56)

matching_sps<-intersect(unique(df.dat.HF$GEN.SPA),unique(db.dat$GEN.SPA)) ## could also add viburnum and vaccinium at the genus level

### we are skipping this to have the most complete data. Models will just completely pool to means on species absant from each other's studies
#df.dat.match<-df.dat.HF %>%filter(GEN.SPA %in% c(matching_sps))  ### make both dataset noly consiston matching sps
#db.dat.match<-db.dat %>%filter(GEN.SPA %in% c(matching_sps))

df.dat.match<-df.dat.HF #These inculde all species
db.dat.match<-db.dat
db.dat.match<-unite(db.dat.match, treatcode, 4:6,sep="")
db.dat.match$Chill<-ifelse(db.dat.match$chilldays==28,0,1)

db.dat.match$FORCE<-ifelse(db.dat.match$temp_day==24,1,0) ### new column where forcing temperature are binary low or high
df.dat.match$FORCE<-ifelse(df.dat.match$warm==20,1,0)


db.dat.match$study<-"DB"
df.dat.match$study<-"DF"

##subset to columns of use
db.dat.match<-dplyr::select(db.dat.match,GEN.SPA,F60,LBB,L11,treatcode,temp_day,photoperiod,Chill,study,FORCE)
df.dat.match<-dplyr::select(df.dat.match,GEN.SPA,fday,bday,lday,treatcode,warm,photo,chill,study,FORCE)




###make columns the same
colnames(db.dat.match)[colnames(db.dat.match)=="F60" ] <- "fday"
colnames(db.dat.match)[colnames(db.dat.match)=="L11" ] <- "lday"
colnames(db.dat.match)[colnames(db.dat.match)=="LBB" ] <- "bday"
colnames(db.dat.match)[colnames(db.dat.match)=="temp_day" ] <- "warm"
colnames(db.dat.match)[colnames(db.dat.match)=="photoperiod" ] <- "photo"
colnames(db.dat.match)[colnames(db.dat.match)=="Chill" ] <- "chill"

colnames(db.dat.match)
colnames(df.dat.match)

unique(db.dat.match$chill)

db.dat.match.short<-filter(db.dat.match,chill==0) ###Use DB's short chill since it is closest in Chilling hours to DF's
###db.dat.match.long<-filter(db.dat.match,chill==1)

both.dat1<-rbind(df.dat.match,db.dat.match) ### since dan F was in he middle of mine, I am going to use both of my chilling
#both.dat2<-rbind(df.dat.match,db.dat.match.long)

both.dat1$PHOTO<-ifelse(both.dat1$photo==8,0,1)
#both.dat2$PHOTO<-ifelse(both.dat2$photo==8,0,1)

both.dat1$warm.cent<-both.dat1$warm-mean(both.dat1$warm)
both.dat1$photo.cent<-both.dat1$photo-mean(both.dat1$photo)
#both.dat2$warm.cent<-both.dat2$warm-mean(both.dat2$warm)
#both.dat2$photo.cent<-both.dat2$photo-mean(both.dat2$photo)


#both.mod.bud1<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat1)
#summary(both.mod.bud1)

#both.mod.bud2<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat2)

#summary(both.mod.bud2)

#both.mod.leaf1<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat1)
#summary(both.mod.leaf1)

#both.mod.leaf2<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(1|GEN.SPA),data=both.dat2)
#summary(both.mod.leaf2)

#both.mod.bud1.int<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat1)
#summary(both.mod.bud1.int)

#both.mod.leaf1.int<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat1)
#summary(both.mod.leaf1.int)

#both.mod.bud2.int<-brm(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat2)
#summary(both.mod.bud2.int)

#both.mod.leaf2.int<-brm(lday~study+FORCE+PHOTO+study:FORCE+study:PHOTO+(bday~study+FORCE+PHOTO+study:FORCE+study:PHOTO|GEN.SPA),data=both.dat2)
#summary(both.mod.leaf2.int)








####but above effects could jsut be cause Dan B used consitantly higher forcing. Luckily my low force average is teh same as DF's high average


#both.dat.sameforce2<-both.dat2 %>% filter(warm %in% c(20,18))# This make a data sheet where the average forcing in both experiments was 15 C (DB low, DF high)
both.dat.sameforce1<-both.dat1 %>% filter(warm %in% c(20,18)) ## Filter just the studies with same forcing

photo.only.bb1<-brm(bday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2500)
summary(photo.only.bb1)
coef(photo.only.bb1)

photo.only.leaf1<-brm(lday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2000)
summary(photo.only.leaf1)

photo.only.flo1<-brm(fday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2000)


fixef(photo.only.bb1, summary = TRUE, robust = FALSE,
     probs = c(0.025, 0.975,.10,.90,.25,.75))

fixef(photo.only.leaf1, summary = TRUE, robust = FALSE,
      probs = c(0.025, 0.975,.10,.90,.25,.75))

fixef(photo.only.flo1, summary = TRUE, robust = FALSE,
      probs = c(0.025, 0.975,.10,.90,.25,.75))
save.image("photovthermoperiod")
#photo.only.bb2<-brm(bday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce2,iter=3000,warmup=2000)
#summary(photo.only.bb2)
#coef(photo.only.bb2)

#photo.only.leaf2<-brm(lday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce2,iter=3000,warmup=2000)
#summary(photo.only.leaf2)



save.image("periodicity_mods.RData")
