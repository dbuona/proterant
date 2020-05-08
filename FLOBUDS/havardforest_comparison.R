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
library(lme4)
library("lmerTest")
library(RColorBrewer)
library(ggstance)
set.seed(613)
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
#HF2<-read.csv("HarvardForest/hf003-06-mean-spp.csv",header=TRUE)
load("fieldexamples")



##Hypothesis 1: variation in FLS is a product of variation in climate between phases
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
weather<-filter(weather,Year>=1990)


weather<-make_all_day_table(weather)

hourtemps<-stack_hourly_temps(weather, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)

979/1440



  HF<-filter(HF,year<=2002)
  HF$treeyear<-paste(HF$tree.id,HF$year)
  HF$burst.flo<-ifelse(!is.na(HF$fopn.jd),HF$fopn.jd,365)
  HF$burst.bb<-ifelse(!is.na(HF$bb.jd),HF$bb.jd,365)

yearios<-unique(HF$treeyear)
df<-data.frame(year=numeric(),GDD.flo=numeric(),treeyear=character())

for(k in c(1:length(yearios))){
gooby<-filter(HF,treeyear==yearios[k])  
goo<-dplyr::filter(hourtemps,Year==gooby$year)  
goo<-filter(goo,JDay<=gooby$burst.flo)
y<-max(GDD(goo$Temp,summ=TRUE,Tbase=5))
dfhere<-data.frame(year=gooby$year,GDD.flo=y,treeyear=gooby$treeyear)
df<-rbind(df,dfhere)
}
HF.flo<-dplyr::select(HF,year,tree.id,burst.flo,treeyear)  
df<-left_join(df,HF.flo)

###now leaves
df2<-data.frame(year=numeric(),GDD.bb=numeric(),treeyear=character())

for(k in c(1:length(yearios))){
  gooby<-filter(HF,treeyear==yearios[k])  
  goo<-dplyr::filter(hourtemps,Year==gooby$year)  
  goo<-filter(goo,JDay<=gooby$burst.bb)
  y<-max(GDD(goo$Temp,summ=TRUE,Tbase=5))
  dfhere2<-data.frame(year=gooby$year,GDD.bb=y,treeyear=gooby$treeyear)
  df2<-rbind(df2,dfhere2)

}

HF.bb<-dplyr::select(HF,year,tree.id,burst.bb,treeyear,species)  
df2<-left_join(df2,HF.bb)

dater<-left_join(df,df2)
dater$GDD.diff<-dater$GDD.bb-dater$GDD.flo
dater$FLS.diff<-dater$burst.bb-dater$burst.flo

###remove the no flos
dater<-filter(dater,burst.flo<365)


write.csv(dater,"GDH_diffs_HF.csv",row.names = FALSE)

#dater<-read.csv(file = "GDH_diffs_HF.csv",header=TRUE)



####simlate field data- heat sums hypothesis

Tb<-5


days<-1:150

rep<-1:20
year<-1989:1999
df<-data.frame(year=factor(),days=numeric(),GDD=numeric(),y=numeric(),tree.id=numeric())

for (k in c(1:length(year))){
  Temp<-rnorm(1,15,5)
    thresh<-rnorm(1,100,3)
    for (j in c(1:length(rep))){
      threshhold<-rnorm(1,thresh,4)
    for (i in c(1:length(days))){
      y<-ifelse((Temp-Tb)*days[i]>=threshhold,days,NA)
      dfhere<-data.frame(year=year[k],days=days[i],GDD=threshhold,y=y,tree.id=rep[j])
      df<-rbind(df,dfhere)
      df<-df[complete.cases(df),]
      
    }}}
df<- df %>%group_by(tree.id,year,GDD) %>%dplyr::summarize(days=min(days))   


scaleFactor <-0.1

a<-ggplot()+geom_point(aes(df$year,df$days),color="#999999",size=.2)+
  stat_summary(aes(df$year,df$days),color="#999999")+
  geom_point(aes(df$year,df$GDD*scaleFactor),color="black",size=.2)+
  stat_summary(aes(df$year,df$GDD*scaleFactor),color="black")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#999999"),
    axis.text.y.left=element_text(color="#999999"),
    axis.title.y.right=element_text(color="black"),
    axis.text.y.right=element_text(color="black"))


days<-1:200

#####now when chilling changes the threshhold for gdd (chilling is latent )
df2<-data.frame(year=factor(),days=numeric(),GDD=numeric(),y=numeric(),tree.id=numeric())

for (k in c(1:length(year))){
  Temp<-rnorm(1,15,5)
  thresh<-rnorm(1,100,30)
  for (j in c(1:length(rep))){
    threshhold<-rnorm(1,thresh,4)
    for (i in c(1:length(days))){
      y<-ifelse((Temp-Tb)*days[i]>=threshhold,days,NA)
      df2here<-data.frame(year=year[k],days=days[i],GDD=threshhold,y=y,tree.id=rep[j])
      df2<-rbind(df2,df2here)
      df2<-df2[complete.cases(df2),]
      
    }}}
df2<- df2 %>%group_by(tree.id,year,GDD) %>%dplyr::summarize(days=min(days))   

#scaleFactor <-mean(df2$days)/ mean(df2$GDD)

b<-ggplot()+geom_point(aes(df2$year,df2$days),color="#999999",size=.2)+
  stat_summary(aes(df2$year,df2$days),color="#999999")+
  geom_point(aes(df2$year,df2$GDD*scaleFactor),color="black",size=.2)+
  stat_summary(aes(df2$year,df2$GDD*scaleFactor),color="black")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#999999"),
    axis.text.y.left=element_text(color="#999999"),
    axis.title.y.right=element_text(color="black"),
    axis.text.y.right=element_text(color="black"))









concept<-ggpubr::ggarrange(a,b)
setEPS()
postscript("..//FLOBUDS/Plots/fieldexamples.eps", width=10, height=4)
concept
dev.off()



###real
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
HF<-filter(HF,year<=2002)
HF$FLS<-HF$bb.jd-HF$fopn.jd
HF$year<-as.factor(HF$year)
Acerub<-filter(dater,species=="ACRU")

scaleFactor <-.1

c<-ggplot(Acerub)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))



jpeg("..//FLOBUDS/Plots/hypothesis1_acerub.jpeg",height=6,width=8,units = "in",res=250)
ggpubr::ggarrange(concept,c,nrow=2,heights=c(.75,1))
dev.off()  


####other species
Acepen<-filter(dater,species=="ACPE")



acpe<-ggplot(Acepen)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))

Qurub<-filter(dater,species=="QURU")


quru<-ggplot(Qurub)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))


Vaco<-filter(dater,species=="VACO")
table(dater$species)

vaco<-ggplot(Vaco)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))

Bepo<-filter(dater,species=="BEPO")


bepo<-ggplot(Bepo)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))


Kala<-filter(dater,species=="KALA")


kala<-ggplot(Kala)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))

Sapu<-filter(dater,species=="SAPU")


sapu<-ggplot(Sapu)+geom_point(aes(as.factor(year),FLS.diff),color="#D55E00",size=1)+
  stat_summary(aes(as.factor(year),FLS.diff),color="#D55E00")+
  geom_point(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73",size=1)+
  stat_summary(aes(as.factor(year),GDD.diff*scaleFactor),color="#009E73")+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./scaleFactor, name="GDDs between phases"))+
  scale_x_discrete("year")+
  ggthemes::theme_base()+
  theme(
    axis.title.y.left=element_text(color="#D55E00"),
    axis.text.y.left=element_text(color="#D55E00"),
    axis.title.y.right=element_text(color="#009E73"),
    axis.text.y.right=element_text(color="#009E73"))

setEPS()
postscript("..//FLOBUDS/Plots/supp_field_sps.eps", width=10, height=11)




ggpubr::ggarrange(acpe,quru,bepo,vaco,kala,sapu,nrow=3,ncol=2,labels = c("a","b","c","d","e","f"))
dev.off()

save.image("fieldexamples")

