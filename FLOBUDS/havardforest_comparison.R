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
weather<-select(weather,c("date","airtmax","airtmin"))
weather<-separate(weather,date,c("Year","Month","Day"),sep="-",remove=TRUE)
colnames(weather)<-c("Year","Month","Day","Tmax","Tmin")
sapply(weather,mode) #mode(weather)
weather$Year<-as.numeric(weather$Year)
weather$Month<-as.numeric(weather$Month)
weather$Day<-as.numeric(weather$Day)

weather<-filter(weather,year>=1989)

all_daylengths<-cbind(JDay=1:365,sapply(daylength(latitude=42.5,JDay=1:365),cbind)) ## calculate day length at HF on every day of year
ad<-as.data.frame(all_daylengths) ## data frame of day length


weather<-make_all_day_table(weather)

hourtemps<-stack_hourly_temps(weather, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)

burst<-Vcraw$bb.jd

##how much cold and GDD did they get before bud burst int aht year
for(k in c(1:length(burst))){ 
  ChillHF<-as.data.frame(chilling(hourtemps,Start_JDay=275,End_JDay=burst[k])) 
}
for(k in c(1:length(burst))){ 
  WarmHF<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=burst[k])) 
}

  
flo<-Vcraw$fopn.jd
for(k in c(1:length(flo))){ 
 ChillHF.flower<-as.data.frame(chilling(hourtemps,Start_JDay=275,End_JDay=flo[k])) 
}

for(k in c(1:length(flo))){ 
  WarmHF.flower<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=flo[k])) 
}


lo<-Vcraw$l75.jd
for(k in c(1:length(lo))){ 
  ChillHF.l75<-as.data.frame(chilling(hourtemps,Start_JDay=275,End_JDay=flo[k])) 
}

for(k in c(1:length(lo))){ 
  WarmHF.l75<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=flo[k])) 
}





WarmHF<-select(WarmHF,c("End_year","Season","GDH")) ## take relevant columns
ChillHF<-select(ChillHF,c("End_year","Season","Chill_portions")) ## take relevant columns

WarmHF.flower<-select(WarmHF.flower,c("End_year","Season","GDH")) ## take relevant columns
ChillHF.flower<-select(ChillHF.flower,c("End_year","Season","Chill_portions"))

WarmHF.l75<-select(WarmHF.l75,c("End_year","Season","GDH")) ## take relevant columns
ChillHF.l75<-select(ChillHF.l75,c("End_year","Season","Chill_portions"))


ChillHF<-filter(ChillHF,End_year %in% c(1990:2001))
WarmHF<-filter(WarmHF,End_year %in% c(1990:2001))

ChillHF.flower<-filter(ChillHF.flower,End_year %in% c(1990:2001))
WarmHF.flower<-filter(WarmHF.flower,End_year %in% c(1990:2001))

ChillHF.l75<-filter(ChillHF.l75,End_year %in% c(1990:2001))
WarmHF.l75<-filter(WarmHF.l75,End_year %in% c(1990:2001))



##make Data sheet
HFtempsbb<-left_join(WarmHF,ChillHF)
HFtempsflo<-left_join(WarmHF.flower,ChillHF.flower)
HFtempsL75<-left_join(WarmHF.l75,ChillHF.l75)
#3make GDD from GDH

HFtempsbb$GDD<-HFtempsbb$GDH/24
HFtempsflo$GDD<-HFtempsflo$GDH/24
HFtempsL75$GDD<-HFtempsL75$GDH/24

colnames(HFtempsbb)<-c("End_year"    ,   "Season"   ,      "GDH.bb"     ,       "Chill_portions.bb" ,"GDD.bb" )
colnames(HFtempsflo)<-c("End_year" ,  "Season","GDH.flo","Chill_portions.flo" ,"GDD.flo")
colnames(HFtempsL75)<-c("End_year" ,  "Season","GDH.L75","Chill_portions.L75" ,"GDD.L75")

HFtempsboth<-left_join(HFtempsbb,HFtempsflo)
HFtempsboth<-left_join(HFtempsboth,HFtempsL75)




#####below might be obsolete

#HFtempsboth$GDDave1<-(HFtempsboth$GDD.bb+HFtempsboth$GDD.flo)/2
#HFtempsboth$CPave1<-(HFtempsboth$Chill_portions.bb+HFtempsboth$Chill_portions.flo)/2

#HFtempsboth$GDDave2<-(HFtempsboth$GDD.L75+HFtempsboth$GDD.flo)/2
#HFtempsboth$CPave2<-(HFtempsboth$Chill_portions.L75+HFtempsboth$Chill_portions.flo)/2



## combine with phenology
mod.dat<-HFtempsboth
colnames(mod.dat)
#select(HFtempsboth,c("End_year","Season","GDDave1","CPave1","GDDave2","CPave2"))
colnames(mod.dat)[1]<-"year"
mod.dat<-left_join(Vcraw,mod.dat,by="year")
##add daylength
colnames(ad)<-c("bb.jd","Sunrise","Sunset","Daylength.leaf")

mod.dat<-left_join(mod.dat,ad,by="bb.jd")
colnames(ad)<-c("fopn.jd","Sunrise","Sunset","Daylength.flo")
mod.dat<-left_join(mod.dat,ad,by="fopn.jd")

colnames(ad)<-c("l75.jd","Sunrise","Sunset","Daylength.l75")
mod.dat<-left_join(mod.dat,ad,by="l75.jd")


#mod.dat$daylenghave1<-(mod.dat$Daylength.leaf+mod.dat$Daylength.flo)/2
#mod.dat$daylenghave2<-(mod.dat$Daylength.l75+mod.dat$Daylength.flo)/2

#mod.daty<-select(mod.dat, c("year"       ,    "tree.id",        "species"  ,      "bb.jd"      ,    "l75.jd"   ,      "fbb.jd"   ,      "fopn.jd",        "FLS"   ,"FLS2"     ,      "Season"      ,   "GDDave1"      ,   "CPave1","GDDave2"      ,   "CPave2","daylenghave1","daylenghave2"))

mod.dat$zchill.bb<-(mod.dat$Chill_portions.bb-mean(mod.dat$Chill_portions.bb,na.rm=TRUE))/sd(mod.dat$Chill_portions.bb,na.rm=TRUE)
ggplot(mod.dat,aes(year,zchill.bb))+geom_bar(stat="identity",position="dodge",fill="blue")

mod.dat$zchill.flo<-(mod.dat$Chill_portions.flo-mean(mod.dat$Chill_portions.flo,na.rm=TRUE))/sd(mod.dat$Chill_portions.flo,na.rm=TRUE)
ggplot(mod.dat,aes(year,zchill.flo))+geom_bar(stat="identity",position="dodge", fill="red")

mod.dat$zphoto.bb<-(mod.dat$Daylength.leaf-mean(mod.dat$Daylength.leaf,na.rm=TRUE))/sd(mod.dat$Daylength.leaf,na.rm=TRUE)
ggplot(mod.dat,aes(year,zphoto.bb))+geom_bar(stat="identity",position="dodge",fill="blue")

mod.dat$zphoto.flo<-(mod.dat$Daylength.flo-mean(mod.dat$Daylength.flo,na.rm=TRUE))/sd(mod.dat$Daylength.flo,na.rm=TRUE)
ggplot(mod.dat,aes(year,zphoto.flo))+geom_bar(stat="identity",position="dodge",fill="red")

mod.dat$zGDD.flo<-(mod.dat$GDD.flo-mean(mod.dat$GDD.flo,na.rm=TRUE))/sd(mod.dat$GDD.flo,na.rm=TRUE)
mod.dat$zGDD.bb<-(mod.dat$GDD.bb-mean(mod.dat$GDD.bb,na.rm=TRUE))/sd(mod.dat$GDD.bb,na.rm=TRUE)

mod.dat$bin_photo_flo<-ifelse(mod.dat$zphoto.flo>0,1,0)
mod.dat$bin_photo_bb<-ifelse(mod.dat$zphoto.bb>0,1,0)

mod.dat$bin_chill_flo<-ifelse(mod.dat$zchill.flo>0,1,0)
mod.dat$bin_chill_bb<-ifelse(mod.dat$zchill.bb>0,1,0)

mod.dat$bin_GDD_flo<-ifelse(mod.dat$zGDD.flo>0,1,0)
mod.dat$bin_GDD_bb<-ifelse(mod.dat$zGDD.bb>0,1,0)



##model
budi<-brm(FLS~bin_chill_bb*bin_GDD_bb,dat=mod.dat)
 summary(budi)
floi<-brm(FLS~bin_chill_flo*bin_GDD_flo,dat=mod.dat)
summary(floi)
###zscore

mod.daty$zchill1<-(mod.daty$CPave1-mean(mod.daty$CPave1,na.rm=TRUE))/sd(mod.daty$CPave1,na.rm=TRUE)
mod.daty$zGDD1<-(mod.daty$GDDave1-mean(mod.daty$GDDave1,na.rm=TRUE))/sd(mod.daty$GDDave1,na.rm=TRUE)
mod.daty$zPhoto1<-(mod.daty$daylenghave1-mean(mod.daty$daylenghave1,na.rm=TRUE))/sd(mod.daty$daylenghave1,na.rm=TRUE)

mod.daty$zchill2<-(mod.daty$CPave2-mean(mod.daty$CPave2,na.rm=TRUE))/sd(mod.daty$CPave2,na.rm=TRUE)
mod.daty$zGDD2<-(mod.daty$GDDave2-mean(mod.daty$GDDave2,na.rm=TRUE))/sd(mod.daty$GDDave2,na.rm=TRUE)
mod.daty$zPhoto2<-(mod.daty$daylenghave2-mean(mod.daty$daylenghave2,na.rm=TRUE))/sd(mod.daty$daylenghave2,na.rm=TRUE)


bbonly<-brm(bb.jd~zGDD1+zchill1+zPhoto1+zGDD1:zchill1+zGDD1:zPhoto1+zchill1:zPhoto1,data=mod.daty)
fonly<-brm(fopn.jd~zGDD1+zchill1+zPhoto1+zGDD1:zchill1+zGDD1:zPhoto1+zchill1:zPhoto1,data=mod.daty)
l75only<-brm(l75.jd~zGDD2+zchill2+zPhoto2+zGDD2:zchill2+zGDD2:zPhoto2+zchill2:zPhoto2,data=mod.daty)


bbonly.nophoto<-brm(bb.jd~zGDD1*zchill1,data=mod.daty)
fonly.nophoto<-brm(fopn.jd~zGDD1*zchill1,data=mod.daty)

extract_cof<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,.25,.75,0.975))),"cue")} 

bbploot<-extract_cof(bbonly)
bbploot$phase<-"leaf budburst"

llploot<-extract_cof(l75only)
llploot$phase<-"leaf expansion"

fploot<-extract_cof(fonly)
fploot$phase<-"flower open"
ploots<-rbind(bbploot,fploot,llploot)
ploots<-dplyr::filter(ploots,cue!="Intercept")

ploots$cue[ploots$cue=="zPhoto2"]<-"zPhoto"
ploots$cue[ploots$cue=="zPhoto1"]<-"zPhoto"
ploots$cue[ploots$cue=="zGDD2"]<-"zGDD"
ploots$cue[ploots$cue=="zGDD1"]<-"zGDD"
ploots$cue[ploots$cue=="zchill1"]<-"zchill"
ploots$cue[ploots$cue=="zchill2"]<-"zchill"
ploots$cue[ploots$cue=="zGDD2:zPhoto2"]<-"zGDD:zPhoto"
ploots$cue[ploots$cue=="zGDD1:zPhoto1"]<-"zGDD:zPhoto"
ploots$cue[ploots$cue=="zchill1:zPhoto1"]<-"zPhoto:zchill"
ploots$cue[ploots$cue=="zchill2:zPhoto2"]<-"zPhoto:zchill"
ploots$cue[ploots$cue=="zGDD1:zchill1"]<-"zGDD:zchill"
ploots$cue[ploots$cue=="zGDD2:zchill2"]<-"zGDD:zchill"
jpeg("..//FLOBUDS/Plots/flo_buds_figures/Vac_at_HF.jpeg",width=1800,height=1000,res=200)
ploots %>%
arrange(Estimate) %>%
  mutate(cue = factor(cue, levels=c("zGDD:zchill",  "zPhoto:zchill","zGDD:zPhoto", "zchill", "zGDD", "zPhoto"))) %>%
ggplot(aes(Estimate,cue))+geom_point(aes(color=phase,shape=phase),size=3)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=phase),size=.5,linetype="dashed",stat="identity",height=0)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),size=.5,height=0)+geom_vline(aes(xintercept=0),color="black")+theme_bw()
dev.off()

ggplot(ploots,aes(Estimate,cue))+geom_point(aes(color=phase,shape=phase),size=3)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=phase),size=.5,linetype="dashed",stat="identity",height=0)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),size=.5,height=0)+geom_vline(aes(xintercept=0),color="black")+theme_bw()
##predict for FLS just
bbtoflo<-brm(FLS~zGDD1+zchill1+zPhoto1+zGDD1:zchill1+zGDD1:zPhoto1+zchill1:zPhoto1,data=mod.daty)

bbtofloplot<-extract_cof(bbtoflo) 



## predictions
##flobuds says flowering is more sensitive to photo
##
bbtofloplot<-filter(bbtofloplot,cue!="Intercept")
jpeg("../FLOBUDS/Plots/flo_buds_figures/effectsize_bbtoflo.jpeg",res=300,width=1500,height=800)
bbtofloplot %>%
  arrange(Estimate) %>%
  mutate(cue = factor(cue, levels=c("zGDD1:zchill1", "zGDD1:zPhoto1", "zchill1:zPhoto1", "zGDD1", "zchill1", "zPhoto1"))) %>%
ggplot(aes(Estimate,cue))+geom_point(size=3)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5),size=.5,linetype="dashed",stat="identity",height=0)+geom_errorbarh(aes(xmin=Q25,xmax=Q75),size=.5,height=0)+geom_vline(aes(xintercept=0),color="black")+ggtitle("Time between budburst and flowering")+theme_bw()
dev.off()
flotolo<-brm(FLS2~zGDD2+zchill2+zPhoto2+zGDD2:zchill2+zGDD2:zPhoto2+zchill2:zPhoto2,data=mod.daty)
flotoloplot<-extract_cof(flotolo)
flotoloplot<-filter(flotoloplot,cue!="Intercept")
jpeg("../FLOBUDS/Plots/flo_buds_figures/effectsize_flotolo.jpeg",res=300,width=1500,height=800)
flotoloplot %>%
  arrange(Estimate) %>%
  mutate(cue = factor(cue, levels=c("zGDD2:zchill2", "zGDD2:zPhoto2", "zchill2:zPhoto2", "zGDD2", "zchill2", "zPhoto2","Intercept"))) %>%
ggplot(aes(Estimate,cue))+geom_point(size=3)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5),size=.5,linetype="dashed",stat="identity",height=0)+geom_errorbarh(aes(xmin=Q25,xmax=Q75),size=.5,height=0)+geom_vline(aes(xintercept=0),color="black")+ggtitle("Time between flowering and L75")
dev.off()

ncol(mod.daty)
plotymody<-gather(mod.daty,measure,accum,17:22)
plotymody$zFLS<-(plotymody$FLS-mean(plotymody$FLS))/sd(plotymody$FLS)
plotymody$zFLS2<-(plotymody$FLS2-mean(plotymody$FLS2))/sd(plotymody$FLS2)

plotymody1<-plotymody %>% filter(measure %in% c("zGDD1","zPhoto1","zchill1"))
plotymody2<-plotymody %>% filter(measure %in% c("zGDD2","zPhoto2","zchill2"))


jpeg("../FLOBUDS/Plots/flo_buds_figures/Havardforestvisual_bbtoflo.jpeg",res=300,width=1500,height=800)
plotymody1  %>%
arrange(year) %>%
mutate(measure = factor(measure, levels=c("zGDD1","zPhoto1","zchill1"))) %>%
ggplot()+theme_bw()+geom_bar(stat="identity",position="dodge",aes(x=year,y=accum,fill=measure))+stat_summary(aes(x=year,y=zFLS),shape=15)+geom_hline(yintercept=0)+ylab("Standardize deviation from mean")
dev.off()

jpeg("../FLOBUDS/Plots/flo_buds_figures/Havardforestvisual_flotolo.jpeg",res=300,width=1500,height=800)
plotymody2  %>%
  arrange(year) %>%
  mutate(measure = factor(measure, levels=c("zGDD2","zPhoto2","zchill2"))) %>%
  ggplot()+theme_bw()+geom_bar(stat="identity",position="dodge",aes(x=year,y=accum,fill=measure))+stat_summary(aes(x=year,y=zFLS2),shape=15)+geom_hline(yintercept=0)+ylab("Standardize deviation from mean")
dev.off()

save.image("HF_comps_vaccor.Rda")
?geom_bar()
