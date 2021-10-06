###compare_periodicity begun by Dan B. on Feb 7 2019
#The goal of this file is to compare the effect of covarying photoperiodicity and thermoperiodcity by comparing dany flynn to Dan b's budburst experients
#### DAn should redo this anaysis
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(brms)
library(chillR)
library(dplyr)
library(tidyr)
library(tibble)
library(ggstance)
library(ggthemes)
library(grid)

setwd("~/Documents/git/proterant/FLOBUDS")
load("periodicity/periodicity_mods.RData")

db.dat<-read.csv("periodicity/Flo_buds_for_thermoperiodicity.csv")
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




both.dat1<-rbind(df.dat.match,db.dat.match) ### since dan F was in he middle of mine, I am going to use both of my chilling
both.dat1$PHOTO<-ifelse(both.dat1$photo==8,0,1)
both.dat1$warm.cent<-both.dat1$warm-mean(both.dat1$warm)
both.dat1$photo.cent<-both.dat1$photo-mean(both.dat1$photo)

##with forcing
goobermody<-brm(bday~study*photo*warm+(study*photo*warm|GEN.SPA),data=both.dat1,iter=3000,warmup=2500)
summary(goobermody)
coef(photo.only.bb1)


both.dat.sameforce1<-both.dat1 %>% filter(warm %in% c(20,18)) ## Filter just the studies with same forcing


photo.only.bb1<-brm(bday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2500)
summary(photo.only.bb1)
coef(photo.only.bb1)
bb.output<-rownames_to_column(as.data.frame(fixef(photo.only.bb1,probs=c(0.1,0.9,0.25,0.75))),"Parameter")

bb.output$Phase<-"budburst"

pd2=position_dodgev(height=0.4)
ggplot(bb.output,aes(Estimate,Parameter))+geom_point(position=pd2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90),linetype="dotted",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="red")+ggtitle("Budburst comparison")

photo.only.leaf1<-brm(lday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce1,iter=3000,warmup=2000)

lo.output<-rownames_to_column(as.data.frame(fixef(photo.only.leaf1,probs=c(0.1,0.9,0.25,0.75))),"Parameter")
lo.output$Phase<-"leafout"
ggplot(lo.output,aes(Estimate,Parameter))+geom_point(position=pd2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90),linetype="dotted",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="red")+ggtitle("leaf out comparison")

veggie<-rbind(bb.output,lo.output)
veggie$Parameter<-ifelse(veggie$Parameter=="studyDF","InterceptDF",veggie$Parameter)
veggie.effect<-veggie %>% filter(Parameter %in% c("PHOTO","studyDF:PHOTO"))

pd2=position_dodgev(height=0.2)
plotyx<-ggplot(veggie,aes(Estimate,Parameter))+geom_point(aes(color=Phase,shape=Phase),position=pd2,size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=Phase),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=Phase),linetype="dashed",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgrey","black"))
plotyy<-ggplot(veggie.effect,aes(Estimate,Parameter))+geom_point(aes(color=Phase, shape=Phase),position=pd2,size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=Phase),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=Phase),linetype="dashed",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgray","black"))+xlim(-10,10)+theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

jpeg("Plots/photothermo_allsps.jpeg",width = 1020, height = 529)
vp <- viewport(width = 0.4, height = 0.5, x = 0.2, y = .99,just=c("left","top"))
plotyx+theme_base()
print(plotyy, vp = vp, )
dev.off()
##3now do matching sps only##########################################################

df.dat.match<-df.dat.HF %>%filter(GEN.SPA %in% c(matching_sps))  ### make both dataset noly consiston matching sps
db.dat.match<-db.dat %>%filter(GEN.SPA %in% c(matching_sps))


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

#db.dat.match.short<-filter(db.dat.match,chill==0) ###Use DB's short chill since it is closest in Chilling hours to DF's
###db.dat.match.long<-filter(db.dat.match,chill==1)


both.dat.match<-rbind(df.dat.match,db.dat.match) ### since dan F was in he middle of mine, I am going to use both of my chilling
#both.dat2<-rbind(df.dat.match,db.dat.match.long)

both.dat.match$PHOTO<-ifelse(both.dat.match$photo==8,0,1)


both.dat.match$warm.cent<-both.dat.match$warm-mean(both.dat.match$warm)
both.dat.match$photo.cent<-both.dat.match$photo-mean(both.dat.match$photo)

both.dat.sameforce.match<-both.dat.match %>% filter(warm %in% c(20,18)) ## Filter just the studies with same forcing

unique(both.dat.sameforce.match$GEN.SPA)

photo.only.bb.match<-brm(bday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce.match,iter=3000,warmup=2500)
photo.only.lo.match<-brm(lday~study*PHOTO+(study*PHOTO|GEN.SPA),data=both.dat.sameforce.match,iter=3000,warmup=2500)

bb.output2<-rownames_to_column(as.data.frame(fixef(photo.only.bb.match,probs=c(0.1,0.9,0.25,0.75))),"Parameter")
bb.output2$Phase<-"budburst"

lo.output2<-rownames_to_column(as.data.frame(fixef(photo.only.lo.match,probs=c(0.1,0.9,0.25,0.75))),"Parameter")
lo.output2$Phase<-"leafout"

veggie2<-rbind(bb.output2,lo.output2)
veggie2$Parameter<-ifelse(veggie2$Parameter=="studyDF","InterceptDF",veggie2$Parameter)

veggie.effect2<-veggie2 %>% filter(Parameter %in% c("PHOTO","studyDF:PHOTO"))


pd2=position_dodgev(height=0.2)
plotyxx<-ggplot(veggie2,aes(Estimate,Parameter))+geom_point(aes(color=Phase,shape=Phase),position=pd2,size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=Phase),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=Phase),linetype="dashed",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgrey","black"))+xlim(-50,80)
plotyyy<-ggplot(veggie.effect2,aes(Estimate,Parameter))+geom_point(aes(color=Phase, shape=Phase),position=pd2,size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=Phase),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=Phase),linetype="dashed",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgray","black"))+xlim(-14,14)+theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

jpeg("Plots/photothermo_matchsps.jpeg",width = 1020, height = 529)
vp <- viewport(width = 0.4, height = 0.5, x = 0.23, y = .99,just=c("left","top"))
plotyxx+theme_base()
print(plotyyy, vp = vp, )
dev.off()

save.image("periodicity_mods.RData")




##a A series of old models
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

