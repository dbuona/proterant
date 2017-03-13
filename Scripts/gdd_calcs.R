## 8 March 2017
# GDD Script and Budburst for HF
# Aim: To find the average accumulated growing degree days for Harvard Forest till budburst
# using John O'Keefe's data

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# Set Working Directory
setwd("~/Documents/git/proterant")
hf<- read.csv("data/WeatherData.csv", header=TRUE)
d<-read.csv("data/hf003-06-mean-spp.csv",header=TRUE)

# Harvard Forest
# To double check my script is accurate
hf<- filter(hf, Site == "hf")
hf$gdd <- hf$AirT - 5
hf$gdd <-ifelse(hf$gdd>0, hf$gdd, 0)
hf$gdd <-ifelse(!is.na(hf$gdd), hf$gdd, 0) #added by dan
hf$count <- ave(
  hf$gdd, hf$Year, 
  FUN=function(x) cumsum(c(0, head(x, -1)))
)
hf<-filter(hf,Year>=1990)

# Clean up dataframes and join
hf<-rename(hf, year=Year)
d<-d%>%
  dplyr::select(year,species, bb.jd,l75.jd,fbb.jd,fopn.jd)
d<- as.data.frame(rapply(object = d, f = round, classes = "numeric", how = "replace", digits = 0)) 
df<-left_join(hf, d)

###species with enough values (from below)
df<-filter(df, species %in% c( "ACPE","ACRU", "ACSA","BEAL","FRAM","QURU"))
##make subset for each phenophase
###budburst
df2<-df%>%
  dplyr::select(year,species, JD, bb.jd, count) %>%
  filter(year>=1990)
df2$day<- ifelse(df2$JD==df2$bb.jd,df2$JD,NA)
df2<-na.omit(df2)
###l75
df3<-df%>%
  dplyr::select(year,species, JD,l75.jd, count) %>%
  filter(year>=1990)
df3$day<- ifelse(df3$JD==df3$l75.jd,df3$JD,NA)
df3<-na.omit(df3)
### fbb.jd
df4<-df%>%
  dplyr::select(year,species, JD, fbb.jd, count) %>%
  filter(year>=1990)
df4$day<- ifelse(df4$JD==df4$fbb.jd,df4$JD,NA)
df4<-na.omit(df4)
##flowers open
df5<-df%>%
  dplyr::select(year,species, JD,fopn.jd, count) %>%
  filter(year>=1990)
df5$day<- ifelse(df5$JD==df5$fopn.jd,df5$JD,NA)
df5<-na.omit(df5)
#plots
#bb
q<-ggplot(df2, aes(x=year, y=count)) + geom_point()
q+facet_wrap(~species)
#l75
q<-ggplot(df3, aes(x=year, y=count)) + geom_point( col="green")
q+facet_wrap(~species)
#fbb
q<-ggplot(df4, aes(x=year, y=count)) + geom_point(col="pink")
q+facet_wrap(~species)
#fopn
q<-ggplot(df5, aes(x=year, y=count)) + geom_point(col="red")
q+facet_wrap(~species)

###is count in each year significantly different? mean and sd 
bb<-df2 %>% group_by(species) %>% summarise(mean(count), sd(count))
l75<-df3 %>% group_by(species) %>% summarise(mean(count), sd(count))
fbb<-df4 %>% group_by(species) %>% summarise(mean(count), sd(count))
fopn<-df5 %>% group_by(species) %>% summarise(mean(count), sd(count))
bb
l75
fbb
fopn
###Things to consider:
#Pattern 1: 1 phenophase considerably more varlaible over time suggest different cues but could also by likely hood to accumulate gdds later in the season
#Pattern 2: if tempurature is main cue, GDD shoudl be minimally varaible. compare varience over time between phenophases within species
#trouble shooting outliers
dfsub<-filter(df,year=="2001")
