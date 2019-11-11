rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")
library(ggplot2)
library(tidyr)
library(dplyr)
library(chillR)
library(brms)
library(tibble)

HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
#HF2<-read.csv("HarvardForest/hf003-06-mean-spp.csv",header=TRUE)
trendovertime.flo<-brm(fopn.jd~year+(year|species),data=HF)  
trendovertime.lo<-brm(l75.jd~year+(year|species),data=HF) 

unique(HF$species)
Vc<-filter(HF,species=="VACO") ## select just vaccinium
#Vc2<-filter(HF2,species=="VACO")
Vcraw<-Vc
Vcraw$FLS<-Vcraw$fopn.jd-Vcraw$bb.jd #time between bb and fopn
Vcraw$FLS2<-Vcraw$l75.jd-Vcraw$fopn.jd # time between fopn and l75
colnames(Vcraw)

### make a plot of phenogy
Vc<-tidyr::gather(Vc, phase,DOY,4:7)
Vc1<-Vc %>% filter(phase %in% c("fopn.jd","l75.jd"))
Vc2<-Vc %>%filter(phase %in% c("fopn.jd","bb.jd"))
ggplot(Vc1,aes(year,DOY))+stat_summary(aes(color=phase))+theme_bw()
ggplot(Vc2,aes(year,DOY))+stat_summary(aes(color=phase))+theme_bw()
## we predict

#1990 91 should be high chilling or forcing, because time between bb and fo is long and fo and l75 is short
#'93 98 and 2000 high chilling or force
#'#94-96 could be low f and c with long photoperoi
#'
#'       

#https://cran.r-project.org/web/packages/chillR/vignettes/hourly_temperatures.html
weather<-read.csv("..//FLOBUDS/data/hf000-01-daily-m.csv",header = TRUE)
weather<-dplyr::select(weather,c("date","airtmax","airtmin"))
weather<-separate(weather,date,c("Year","Month","Day"),sep="-",remove=TRUE)
colnames(weather)<-c("Year","Month","Day","Tmax","Tmin")
sapply(weather,mode) #mode(weather)
weather$Year<-as.numeric(weather$Year)
weather$Month<-as.numeric(weather$Month)
weather$Day<-as.numeric(weather$Day)

unique(weather$Year)
weather<-filter(weather,Year>=1989)

all_daylengths<-cbind(JDay=1:365,sapply(daylength(latitude=42.5,JDay=1:365),cbind)) ## calculate day length at HF on every day of year
ad<-as.data.frame(all_daylengths) ## data frame of day length


weather<-make_all_day_table(weather)

hourtemps<-stack_hourly_temps(weather, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)

burst<-Vcraw$bb.jd

##how much cold and GDD did they get before bud burst int aht year
for(k in c(1:length(burst))){ 
  ChillHF<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=60)) 
}
for(k in c(1:length(burst))){ 
  WarmHF<-as.data.frame(chilling(hourtemps,Start_JDay=60,End_JDay=burst[k])) 
}

  
flo<-Vcraw$fopn.jd
for(k in c(1:length(flo))){ 
 ChillHF.flower<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=60)) 
}

for(k in c(1:length(flo))){ 
  WarmHF.flower<-as.data.frame(chilling(hourtemps,Start_JDay=60,End_JDay=flo[k])) 
}


#lo<-Vcraw$l75.jd
#for(k in c(1:length(lo))){ 
 # ChillHF.l75<-as.data.frame(chilling(hourtemps,Start_JDay=275,End_JDay=lo[k])) 
#}

#for(k in c(1:length(lo))){ 
 # WarmHF.l75<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=lo[k])) 
#}





WarmHF<-dplyr::select(WarmHF,c("End_year","Season","GDH")) ## take relevant columns
ChillHF<-dplyr::select(ChillHF,c("End_year","Season","Chill_portions")) ## take relevant columns

WarmHF.flower<-dplyr::select(WarmHF.flower,c("End_year","Season","GDH")) ## take relevant columns
ChillHF.flower<-dplyr::select(ChillHF.flower,c("End_year","Season","Chill_portions"))

#WarmHF.l75<-dplyr::select(WarmHF.l75,c("End_year","Season","GDH")) ## take relevant columns
#ChillHF.l75<-dplyr::select(ChillHF.l75,c("End_year","Season","Chill_portions"))


ChillHF<-filter(ChillHF,End_year %in% c(1990:2001))
WarmHF<-filter(WarmHF,End_year %in% c(1990:2001))
ChillHF.flower<-filter(ChillHF.flower,End_year %in% c(1990:2001))
WarmHF.flower<-filter(WarmHF.flower,End_year %in% c(1990:2001))

#ChillHF.l75<-filter(ChillHF.l75,End_year %in% c(1990:2001))
#WarmHF.l75<-filter(WarmHF.l75,End_year %in% c(1990:2001))



##make Data sheet
HFtempsbb<-left_join(WarmHF,ChillHF)
HFtempsflo<-left_join(WarmHF.flower,ChillHF.flower)
#HFtempsL75<-left_join(WarmHF.l75,ChillHF.l75)
#3make GDD from GDH

HFtempsbb$GDD<-HFtempsbb$GDH/24
HFtempsflo$GDD<-HFtempsflo$GDH/24
#HFtempsL75$GDD<-HFtempsL75$GDH/24

colnames(HFtempsbb)
HFtempsbb$GDD.z<-(HFtempsbb$GDD-mean(HFtempsbb$GDD))/sd(HFtempsbb$GDD)
HFtempsbb$Chill_portions.z<-(HFtempsbb$Chill_portions-mean(HFtempsbb$Chill_portions,na.rm=TRUE))/sd(HFtempsbb$Chill_portions,na.rm=TRUE)

HFtempsbb.plot<-gather(HFtempsbb,factor,value,6:7)
##1990,1992,1993,1999
ggplot(HFtempsbb.plot,aes(End_year,value,fill=factor))+geom_bar(stat="identity",position="dodge")+theme_linedraw()
#1992,1993,1994,1997,1999

#check if its the same years with flowering
HFtempsflo$GDD.z<-(HFtempsflo$GDD-mean(HFtempsflo$GDD))/sd(HFtempsflo$GDD)
HFtempsflo$Chill_portions.z<-(HFtempsflo$Chill_portions-mean(HFtempsflo$Chill_portions,na.rm=TRUE))/sd(HFtempsflo$Chill_portions,na.rm=TRUE)
HFtempsflo.plot<-gather(HFtempsflo,factor,value,6:7)
ggplot(HFtempsflo.plot,aes(End_year,value,fill=factor))+geom_bar(stat="identity",position="dodge")+theme_linedraw()
#1992,1993,1997,1999


colnames(HFtempsbb)<-c("End_year"    ,   "Season"   ,"GDH.bb",  "CP.bb", "GDD.bb","GDD.bb.z","CP.z" )
colnames(HFtempsflo)<-c("End_year"    ,   "Season"   ,"GDH.flo",  "CP.flo", "GDD.flo","GDD.flo.z","CP.z" )
#colnames(HFtempsL75)<-c("End_year" ,  "Season","GDH.L75","Chill_portions.L75" ,"GDD.L75")

###



HFtempsboth<-left_join(HFtempsbb,HFtempsflo)
HFtempsboth<-filter(HFtempsboth,End_year %in% c(1992,1993,1994,1997,1999))
bb<-brm(GDD.bb.z~CP.z,data=HFtempsboth)
flo<-brm(GDD.flo.z~CP.z,data=HFtempsboth)

summary(bb)
summary(flo)
#HFtempsboth<-left_join(HFtempsboth,HFtempsL75)









alldats %>%
  arrange(Estimate) %>%
  mutate(cue = factor(cue, levels=c("chill portions"))) %>%
  ggplot(aes(Estimate,cue))+geom_point(size=3,aes(color=phase))+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=phase),size=.5,linetype="dashed",stat="identity",height=0)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),size=.5,height=0)+geom_vline(aes(xintercept=0),color="black")+ggtitle("Time between flowering and L75")




save.image("HF_comps_vaccor.Rda")
?geom_bar()
