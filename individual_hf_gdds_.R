###for indv species
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(lme4)
library(car)
library(arm)
#install.packages("broom")
library(broom)
library(tidyr)
setwd("~/Documents/git/proterant")
hf<- read.csv("data/WeatherData.csv", header=TRUE)
hf<-distinct(hf,Date,.keep_all = TRUE)
d<-read.csv("data/hf003-05-mean-ind.csv",header=TRUE)
###adding gdd
hf$gdd <- hf$AirT - 5
hf$gdd <-ifelse(hf$gdd>0, hf$gdd, 0)
hf$gdd <-ifelse(!is.na(hf$gdd), hf$gdd, 0) #added by dan-this is the problem line
hf$count <- ave(
  hf$gdd, hf$Year, 
  FUN=function(x) cumsum(c(0, head(x, -1)))
)
hf<-filter(hf,Year>=1990)
###mergine datasets
hf<-rename(hf, year=Year)
#d<-d%>%
  #dplyr::select(year,species, bb.jd,l75.jd,fbb.jd,fopn.jd)
d<- as.data.frame(rapply(object = d, f = round, classes = "numeric", how = "replace", digits = 0)) 
df<-left_join(hf, d)

##subsetting floral and leaf bud burst
###budburst
df2<-df%>%
  dplyr::select(year,species, JD, bb.jd, tree.id,count)
df2$day<- ifelse(df2$JD==df2$bb.jd,df2$JD,NA)
df2<-na.omit(df2)
q<-ggplot(df2, aes(x=year, y=count)) + geom_point()+ggtitle("bud burst")
q+facet_wrap(~species)
###fbb
df3<-df%>%
  dplyr::select(year,species, JD, fbb.jd, tree.id,count)
df3$day<- ifelse(df3$JD==df3$fbb.jd,df3$JD,NA)
df3<-na.omit(df3)
q<-ggplot(df3, aes(x=year, y=count)) + geom_point(aes(col="pink"))+ggtitle("floral bud burst")
q+facet_wrap(~species)

recap<-df2 %>%
  group_by(tree.id,species) %>%
  summarise(mean.bb = mean(count),sd.bb=sd(count))
#now divid sd, by mean

recap<-mutate(recap, adj.bb= sd.bb/mean.bb)
recapf<-df3 %>%
  group_by(tree.id,species) %>%
  summarise(mean.fbb = mean(count),sd.fbb=sd(count))
recapf<-mutate(recapf, adj.fbb= sd.fbb/mean.fbb)
budbu<-left_join(recap, recapf)

#budbu<-gather(budbu,class, gddtoevent, sd.bb:sd.fbb)
#list2env(split(budbu, budbu$species), envir = .GlobalEnv)
##are bb and fbb sd's correlated?
list_sp<- split(budbu, budbu$species)
sapply(list_sp, function(x) cor(x$adj.bb,x$adj.fbb))
t.test(budbu$sd.bb,budbu$sd.fbb)
###need to figure out a way to do this for each species... also remind me why i needed to standardize the sd with adj?
