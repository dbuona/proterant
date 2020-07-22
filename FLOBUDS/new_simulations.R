rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS/")

library(eesim)
library(lubridate)
library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(broom)

#simulate a year of data

weather<-read.csv("data/hf000-01-daily-m.csv")
weather<-dplyr::select(weather,date,airt)
weather$heatsum<-weather$airt-5
weather$heatsum<-ifelse(weather$heatsum>=0,weather$heatsum,0)
weather$heatsum<-ifelse(is.na(weather$heatsum),0,weather$heatsum)
weather$doy<-yday(weather$date)
weather<-tidyr::separate(weather, date, c("year","month","day"),sep="-")


weather<-weather %>% group_by(year)%>% mutate(GDD = cumsum(heatsum))

weather<-filter(weather,year>=1990)

year<-1990:2002
inds<-1:20
#df<-data.frame(year=character(),flowering=numeric(),leafing=numeric())

#for(i in c(1:length(year))){
#dof<-rnorm(1,50,1)
#dol<-rnorm(1,80,1)
#dfhere<-data.frame(year=as.character(year[i]),GDD.flo=dof,GDD.leaf=dol)

#df<-rbind(df,dfhere)}
#weather<-left_join(weather,df)
  

df<-data.frame(year=character(),tree.id=numeric(),chill=numeric(),flowering=numeric(),leafing=numeric())

for(i in c(1:length(year))){
  chill<-sample(c(0,1),1)
for(j in c(1:length(inds))){
  dof<-rnorm(1,40,1)-(rnorm(1,3,.5)*chill)
  dol<-rnorm(1,80,1)-(rnorm(1,3,.5)*chill)
  dfhere<-data.frame(year=as.character(year[i]),tree.id=inds[j],chill=chill,GDD.flo=dof,GDD.leaf=dol)
  
  df<-rbind(df,dfhere)}}
weather<-left_join(weather,df)

weather.flow<-weather%>% group_by(year,tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) 
weather.leaf<-weather%>% group_by(year,tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) 



weather<-cbind(weather.leaf,weather.flow)
weather$diff.doy<-weather$doy-weather$doy1
weather$diff.GDD<-weather$GDD.leaf-weather$GDD.flo

a<-ggplot()+stat_summary(data=weather,aes(year,diff.doy),color="black")+
  stat_summary(data=weather,aes(year,diff.GDD*1),color="gray")+ggthemes::theme_base(base_size = 10)+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./1, name="GDDs between phases"))+
  theme(
    axis.title.y.right=element_text(color="#999999"),
    axis.text.y.right=element_text(color="#999999")
  )


#simulate a year of data

weather2<-read.csv("data/hf000-01-daily-m.csv")
weather2<-dplyr::select(weather2,date,airt)
weather2$heatsum<-weather2$airt-5
weather2$heatsum<-ifelse(weather2$heatsum>=0,weather2$heatsum,0)
weather2$heatsum<-ifelse(is.na(weather2$heatsum),0,weather2$heatsum)
weather2$doy<-yday(weather2$date)
weather2<-tidyr::separate(weather2, date, c("year","month","day"),sep="-")


weather2<-weather2 %>% group_by(year)%>% mutate(GDD = cumsum(heatsum))

weather2<-filter(weather2,year>=1990)

year<-1990:2002
inds<-1:20
#df<-data.frame(year=character(),flowering=numeric(),leafing=numeric())

#for(i in c(1:length(year))){
#dof<-rnorm(1,50,1)
#dol<-rnorm(1,80,1)
#dfhere<-data.frame(year=as.character(year[i]),GDD.flo=dof,GDD.leaf=dol)

#df<-rbind(df,dfhere)}
#weather2<-left_join(weather2,df)




####3 

df<-data.frame(year=character(),tree.id=numeric(),chill=numeric(),flowering=numeric(),leafing=numeric())

for(i in c(1:length(year))){
  chill<-sample(c(0,1,2,3,4,5),1)
  for(j in c(1:length(inds))){
    dof<-rnorm(1,40,1)-(rnorm(1,3,.5)*chill)
    dol<-rnorm(1,80,1)-(rnorm(1,6,.5)*chill)
    dfhere<-data.frame(year=as.character(year[i]),tree.id=inds[j],chill=chill,GDD.flo=dof,GDD.leaf=dol)
    
    df<-rbind(df,dfhere)}}
weather2<-left_join(weather2,df)

weather2.flow<-weather2%>% group_by(year,tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) 
weather2.leaf<-weather2%>% group_by(year,tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) 



weather2<-cbind(weather2.leaf,weather2.flow)
weather2$diff.doy<-weather2$doy-weather2$doy1
weather2$diff.GDD<-weather2$GDD.leaf-weather2$GDD.flo

b<-ggplot()+stat_summary(data=weather2,aes(year,diff.doy),color="black")+
  stat_summary(data=weather2,aes(year,diff.GDD*1),color="gray")+ggthemes::theme_base(base_size = 10)+
  scale_y_continuous("Days between phases", sec.axis=sec_axis(~./1, name="GDDs between phases"))+
  theme(
    axis.title.y.right=element_text(color="#999999"),
    axis.text.y.right=element_text(color="#999999")
  )

ggpubr::ggarrange(a,b)






