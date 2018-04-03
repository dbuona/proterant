###This is the main cleaning file for flobuds data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS")
d<-read.csv("datasheet_alt_format.csv",header=TRUE)
library(tidyverse)
library(lubridate)
library("Hmisc")

###Clean date
unique(d$date)
d$date[d$date=="26-Nov"]<-"11/26/17"
d$Date<-d$date
d<-separate(d,date,c("month","day","year"))

d$Date<-as.Date(d$Date,format =  "%m/%d/%y")

d$doy<-yday(d$Date)
unique(d$doy)
start<-yday("2017/11/21")
365-325
d$doy.adjusted<-ifelse(d$year==17,d$doy-start,40+(d$doy))
unique(d$doy.adjusted)

#cold treatment wen in Dec/19/2017
start2<-yday("2017/12/19")
365-353
d$doy.adjusted2<-ifelse(d$year==17,d$doy-start2,12+(d$doy))
unique(d$doy.adjusted2)

d$doy.final<-ifelse(d$Chill==0,d$doy.final<-d$doy.adjusted,d$doy.final<-d$doy.adjusted2)

d<-filter(d,doy.final>=0)
unique(d$doy.final)
max(d$doy.final)
q<-filter(d,Chill==1)
r<-filter(d,Chill==0)
unique(q$doy.final)
unique(r$doy.final)
###give each entry a unique id
d<-unite(d,id,name,good_flaskid,sep="_",remove = FALSE)

###clean treatment
d<-unite(d,treatcode,Force,Light,Chill,sep="_",remove = FALSE)
d$Light<-ifelse(d$treatcode=="W_S_0","L",d$Light)
d$Light<-ifelse(d$treatcode=="W_L_0","S",d$Light)
d<-dplyr::select(d,-treatcode)## we'll add this later, but in the meantime it will mess up all the gather() commands if we dont drop it
#### gather locations so it is informative but not catagorical
dx<-gather(d,"flower_location","flophase",c(8,10))
dx<-gather(dx,"leaf_location","leafphase",c(8,9))

####clean flowering
unique(dx$flophase)####this needs alot of cleaning!!!!
###79 flowering catagories to clean
###change all to format: mixed,female,male
dx$flophase[dx$flophase=="0"]<-"0,0,0"
dx$flophase[dx$flophase=="-"]<-"0,0,0"
dx$flophase[dx$flophase=="9"]<-"9,0,0"
dx$flophase[dx$flophase=="7"]<-"7,0,0"
dx$flophase[dx$flophase=="1"]<-"1,0,0"
dx$flophase[dx$flophase=="NA"]<-"NA,NA,NA"
dx$flophase[is.na(dx$flophase)]<-"NA,NA,NA"
dx$flophase[dx$flophase==""]<-"0,0,0"
dx$flophase[dx$flophase=="na"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="n"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="55"]<-"55,0,0"
dx$flophase[dx$flophase=="51"]<-"51,0,0"
dx$flophase[dx$flophase=="60"]<-"60,0,0"
dx$flophase[dx$flophase=="61"]<-"61,0,0"
dx$flophase[dx$flophase=="63"]<-"63,0,0"
dx$flophase[dx$flophase=="64"]<-"64,0,0"
dx$flophase[dx$flophase=="62"]<-"62,0,0"
dx$flophase[dx$flophase=="65"]<-"65,0,0"
dx$flophase[dx$flophase=="60-M"]<-"0,0,60-M"
dx$flophase[dx$flophase=="67"]<-"0,0,67"
dx$flophase[dx$flophase=="1-M"]<-"0,0,1-M"
dx$flophase[dx$flophase=="51-M"]<-"0,0,51-M"
dx$flophase[dx$flophase=="X"]<-"X,X,X"
dx$flophase[dx$flophase=="55-M"]<-"0,0,55-M"
dx$flophase[dx$flophase=="67-M"]<-"0,0,67-M"
dx$flophase[dx$flophase=="69"]<-"69,0,0"
dx$flophase[dx$flophase=="62-M"]<-"0,0,62-M"
dx$flophase[dx$flophase=="65-M"]<-"0,0,65-M"
dx$flophase[dx$flophase=="65-M,60-F"]<-"0,60-F,65-M"
dx$flophase[dx$flophase=="67,9"]<-"67,0,0"
dx$flophase[dx$flophase=="69,51"]<-"69,0,0"
dx$flophase[dx$flophase=="65-M,65-F"]<-"0,65-F,65-M"
dx$flophase[dx$flophase=="69,9"]<-"69,0,0"
dx$flophase[dx$flophase=="69,55"]<-"69,0,0"
dx$flophase[dx$flophase=="67-M,F"]<-"0,67-F,67-M"
dx$flophase[dx$flophase=="69-M"]<-"0,0,69-M"
dx$flophase[dx$flophase=="667-M"]<-"0,0,67-M"
dx$flophase[dx$flophase=="11"]<-"11,0,0"
dx$flophase[dx$flophase=="D"]<-"D,D,D"
dx$flophase[dx$flophase=="63-M"]<-"0,0,63-M"
dx$flophase[dx$flophase=="60-F"]<-"0,60-F,0"
dx$flophase[dx$flophase=="60-F"]<-"0,60-F,0"
dx$flophase[dx$flophase=="Na"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="69-"]<-"69?,69?,?"
dx$flophase[dx$flophase=="lost"]<-"0,0,lost"
dx$flophase[dx$flophase=="60-m"]<-"0,0,60-M"
dx$flophase[dx$flophase=="?"]<-"?,?,?"
dx$flophase[dx$flophase=="61-M"]<-"0,0,61-M"
dx$flophase[dx$flophase=="60-F"]<-"0,59-F,0"
dx$flophase[dx$flophase=="55-F"]<-"0,55-F,0"
dx$flophase[dx$flophase=="61-F"]<-"0,61-F,0"
dx$flophase[dx$flophase=="63-F"]<-"0,63-F,0"
dx$flophase[dx$flophase=="65-F"]<-"0,65-F,0"
dx$flophase[dx$flophase=="65%"]<-"65,0,0"
dx$flophase[dx$flophase=="51-F"]<-"0,51-F,0"
dx$flophase[dx$flophase=="62-F"]<-"0,62-F,0"
dx$flophase[dx$flophase=="51?"]<-"51?,0,0"
dx$flophase[dx$flophase=="67-F"]<-"0,67-F,0"
dx$flophase[dx$flophase=="65-F,55-M"]<-"0,65-F,55-M"
dx$flophase[dx$flophase=="N"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="67-F,63-M"]<-"0,67-F,63-M"
dx$flophase[dx$flophase=="67-F,65-M"]<-"0,67-F,65-M"
dx$flophase[dx$flophase=="69-F"]<-"0,69-F,0"
dx$flophase[dx$flophase=="69-F,51-M"]<-"0,69-F,51-M"
dx$flophase[dx$flophase=="69-F,67-M"]<-"0,69-F,67-M"
dx$flophase[dx$flophase=="65-F,51-M"]<-"0,65-F,51-M"
dx$flophase[dx$flophase=="69-F,51"]<-"0,69-F,51-M"
dx$flophase[dx$flophase=="69-F,67"]<-"0,69-F,67-M"
dx$flophase[dx$flophase=="69-M,60-F"]<-"0,60-F,69-M"
dx$flophase[dx$flophase=="65-m"]<-"0,0,65-M"
dx$flophase[dx$flophase=="69-F,60-M"]<-"0,69-F,60-M"
dx$flophase[dx$flophase=="69-F,69-M"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="60-M,69-F"]<-"0,69-F,60-M"
dx$flophase[dx$flophase=="69-M,65-F"]<-"0,65-F,69-M"
dx$flophase[dx$flophase=="67-F,60-M"]<-"0,67-F,60-M"
dx$flophase[dx$flophase=="60-F,55-M"]<-"0,60-F,55-M"
dx$flophase[dx$flophase=="69-M,69-F"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="69-M,67-F"]<-"0,67-F,69-M"
dx$flophase[dx$flophase=="69-M,57-M"]<-"0,67-F,69-M"
dx$flophase[dx$flophase=="69-F,65-M"]<-"0,69-F,65-M"
dx$flophase[dx$flophase=="69-F,65-M"]<-"0,69-F,65-M"
dx$flophase[dx$flophase=="15"]<-"15?,0,0"
dx$flophase[dx$flophase=="59-F"]<-"0,59-F,0"
unique(dx$flophase)
###more cleaning
dx$flophase[dx$flophase=="0,69-F,67-F"]<-"0,69-F,67-M"
dx$flophase[dx$flophase=="0,065-F,0"]<-"0,65-F,0"
dx$flophase[dx$flophase=="0,60-F,65-N"]<-"0,60-F,65-M"
dx$flophase[dx$flophase=="0,60-F,060-M"]<-"0,60-F,60-M"
dx$flophase[dx$flophase=="0,69-F,690M"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="0,69-M,69-M"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="0,60-F,60M"]<-"0,60-F,60-M"
dx$flophase[dx$flophase=="0,69-F,065-M"]<-"0,69-F,65-M"
unique(dx$flophase)


dx<-separate(dx,flophase,c("mixphase","femphase","malephase"),sep=",")
dx<-gather(dx,flotype,flophase,17:19)
unique(dx$flophase)
###### make everything 60######## for computation sake this way if first flower was score at 65 (etc) it makes the analysis.
dx$flophase[dx$flophase=="61"]<-"60"
dx$flophase[dx$flophase=="62"]<-"60"
dx$flophase[dx$flophase=="63"]<-"60"
dx$flophase[dx$flophase=="64"]<-"60"
dx$flophase[dx$flophase=="65"]<-"60"
dx$flophase[dx$flophase=="67"]<-"60"

dx$flophase[dx$flophase=="61-F"]<-"60-F"
dx$flophase[dx$flophase=="62-F"]<-"60-F"
dx$flophase[dx$flophase=="63-F"]<-"60-F"
dx$flophase[dx$flophase=="64-F"]<-"60-F"
dx$flophase[dx$flophase=="65-F"]<-"60-F"
dx$flophase[dx$flophase=="67-F"]<-"60-F"

dx$flophase[dx$flophase=="61-M"]<-"60-M"
dx$flophase[dx$flophase=="62-M"]<-"60-M"
dx$flophase[dx$flophase=="63-M"]<-"60-M"
dx$flophase[dx$flophase=="64-M"]<-"60-M"
dx$flophase[dx$flophase=="65-M"]<-"60-M"
dx$flophase[dx$flophase=="67-M"]<-"60-M"
dx$flophase[dx$flophase=="60M"]<-"60-M"
table(dx$flophase)

### clean flophase
dx$flotype<-ifelse(dx$GEN.SPA=="VAC.COR","mixphase",dx$flotype)
dx$flotype<-ifelse(dx$GEN.SPA=="PRUN.PEN","mixphase",dx$flotype)
dx$flotype<-ifelse(dx$GEN.SPA=="PRU.VIR","mixphase",dx$flotype)

#to do:
#clean individuals amalanchier and bet lenta snuck in. Make sure that is accounted for
#clean the few COM.PERs that are reported as mixphase
#put in treatment values (IE 8 and 12 photoperiod)
#Calculate chilling

#############CREATE A DATA SHEET THAT HAS THE THREE FLOWER TYPES IN SEPARATE COLUMNS#########################
#################THIS WILL BE USEFUL FOR COMPARING CHANGES IN PROTANDRY OR PROTOGYNY######################

####find the first day when species reached 15
d.leaf<-filter(dx,leafphase==15)
first<-aggregate(d.leaf$doy.final, by = list(d.leaf$id), min)
####combine with all data
dater<-as.data.frame(unique(dx$id))
colnames(dater)<- c("id")
colnames(first)<-c("id","leaf_day")
dater<-full_join(dater,first,by="id") ###now you have a data set with first leaves
### flowers (currrently mixed only)
d.flo<-filter(dx,flophase==60)
firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id), min)
colnames(firstflo)<-c("id","flo_day")
dater<-full_join(dater,firstflo,by="id")
###add female
d.flo<-filter(dx,flophase=="60-F")
firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id), min)
colnames(firstflo)<-c("id","flo_dayF")
dater<-full_join(dater,firstflo,by="id")
###add male
d.flo<-filter(dx,flophase=="60-M")
firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id), min)
colnames(firstflo)<-c("id","flo_dayM")
dater<-full_join(dater,firstflo,by="id")
##### pull experimenta; information to merge. THis give a full data sheet 
###in which male, female and mixd flowering can be compared
dxx<-select(dx,1:7)
dxx<-distinct(dxx)
good.dat<-left_join(dater,dxx)
good.dat<-gather(good.dat,sex,DOY,2:5)
good.dat<-unite(good.dat,treatment,Force,Light,Chill,sep="" )

#ggplot(good.dat,aes(x=treatment,y=flo_day_Mon, color=sex))+geom_point()+facet_wrap(~GEN.SPA)
ggplot(good.dat,aes(x=treatment,y=DOY))+stat_summary(aes(color=sex))+geom_point(size=0.5,aes(color=sex))+facet_wrap(~GEN.SPA)
##########################################################################################

############MAKE A DATE SHEET THAT AGGREGATES ALL FLOWER TO EVALUATE COARSE HYSTERANTHY############
DAT<-separate(dx,flophase,c("abosolute_flower","sex_ind"),sep="-")
L<-filter(DAT,leafphase==15)
L1<-aggregate(L$doy.final, by = list(L$id), min)

datUM<-as.data.frame(unique(DAT$id))
colnames(datUM)<- c("id")
colnames(L1)<-c("id","leaf_day")
datUM<-full_join(datUM,L1,by="id") ###now you have a data set with first leaves

### flowers (currrently mixed only)
Fl<-filter(DAT,abosolute_flower==60)
Fl1<-aggregate(Fl$doy.final, by = list(Fl$id), min)
colnames(Fl1)<-c("id","flo_day")
datUM<-full_join(datUM,Fl1,by="id")
great.dat<-left_join(datUM,dxx)
great.dat<-unite(great.dat,treatment,Force,Light,Chill,sep="" )
#####################################
###good.dat, great.dat are the files to analyze
########################################################################
### metrics for great date##############################################
###how many entries have full entries?
full<-subset(great.dat, !is.na(great.dat$leaf_day)&!is.na(great.dat$flo_day))
nrow(full)#128 our of 576
table(full$GEN.SPA)
table(full$treatment)
### how many have flowering?
flowerfun<-subset(great.dat,!is.na(great.dat$flo_day))
nrow(flowerfun) #141 out of 576
table(flowerfun$GEN.SPA)
table(flowerfun$treatment)

#one or the other:
something<-subset(great.dat, !is.na(great.dat$leaf_day)|!is.na(great.dat$flo_day))
nrow(something) #352
352/576
table(something$GEN.SPA)
table(something$treatment)

############################################################################################################
###plot the raw data

an.data<-gather(great.dat,Phenophase,DOY,2:3)


bigsp<-filter(an.data, GEN.SPA %in% c( "COM.PER","COR.COR","ILE.MUC", "PRU.PEN","VAC.COR"))
ggplot(bigsp, aes(x=treatment, y=DOY,color=Phenophase))+stat_summary()+geom_point(size=.25)+facet_wrap(~GEN.SPA)

berries<-filter(an.data, GEN.SPA=="VAC.COR")
ggplot(berries, aes(x=treatment, y=DOY,color=Phenophase))+stat_summary()+geom_point(size=.25)
ggplot(an.data, aes(x=treatment, y=DOY,))+geom_point(aes(color=Phenophase))+facet_wrap(~GEN.SPA)
