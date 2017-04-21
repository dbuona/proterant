##Disperate bud behaviors floral and leaf phenophased explored using long term dataset from Harvard Forest (O'Keefe).
##This script compiles and organizes several older scripts which on completion will be move to folder "old.scripts"
#Began 3/30/2017

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
library(plyr)
#install.packages("broom")
library(broom)
library(tidyr)
setwd("~/Documents/git/proterant")
hf<- read.csv("data/WeatherData.csv", header=TRUE)
hf<-distinct(hf,Date,.keep_all = TRUE)
d<-read.csv("data/hf003-05-mean-ind.csv",header=TRUE)
d.spp<-read.csv("data/hf003-06-mean-spp.csv",header = TRUE)
##################################################################
#I. Compare trends in flowering and leafout timing over the duration of the study 1990-2011
###currently using means calculated by O'Keefe (datasheet d.spp)
##plotting
#d.1<-filter(d.spp, species %in% c( "ACPE","ACRU", "ACSA","BEAL","FRAM","QURU")) #filter species where observation go full legnth of study
d.1<-d.spp
d.2<-gather(d.1,phenophase, eventday, bb.jd:fopn.jd)

#d.1<-d.1 %>% group_by(year, species,phenophase) %>% summarise(meaneventday=mean(eventday))
##plotting all species
ggplot(d.2, aes(year, eventday, colour = factor(phenophase))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
###slopes for species
slopes.spp<-d.2 %>% 
  group_by(phenophase) %>% 
  do({
    mod = lm(eventday ~ year, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

##now plotting individual species
ggplot(d.2, aes(year, eventday, colour = factor(phenophase))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)+facet_wrap(~species)
##now slopes
slopes.ind<-d.2 %>% 
  group_by(phenophase,species) %>% 
  do({
    mod = lm(eventday ~ year, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
slopes.spp
slopes.ind
##################Great#############################################################
#Now let's compare the how the offset between leaf and bud proccesses change year to year
##########################################################################################
##calculate annual offset between fbb-bb, fopn-l75 and bb-fopn
offset<-mutate(d.1, buds_seq = bb.jd -fbb.jd)
offset<-mutate(offset, functional_seq = l75.jd -fopn.jd)
offset<-mutate(offset, early_seq = bb.jd -fopn.jd)

offset<-group_by(offset, year,species)
###plot flower vs. leaf offsets over time
ggplot(offset, aes(year,buds_seq))+ geom_point()+(geom_abline(intercept=0, slope=0))+facet_wrap(~species)+ggtitle("Leaf vs. Flower Budbust")+geom_smooth(method = "lm", se = FALSE)
ggplot(offset, aes(year,functional_seq))+ geom_point()+(geom_abline(intercept=0, slope=0))+facet_wrap(~species)+ggtitle("Leaf vs. Flower Developement")+geom_smooth(method = "lm", se = FALSE)
ggplot(offset, aes(year,early_seq))+ geom_point()+(geom_abline(intercept=0, slope=0))+facet_wrap(~species)+ggtitle("Leaf Budburst vs. Flower Developement")+geom_smooth(method = "lm", se = FALSE)
#####Would this look different with multiple indivudals rather than species means data
d.3<-filter(d, species %in% c( "ACPE","ACRU", "ACSA","BEAL","FRAM","QURU"))
offset2<-mutate(d.3, buds_seq = bb.jd -fbb.jd)
offset2<-mutate(offset2, functional_seq = l75.jd -fopn.jd)
offset2<-mutate(offset2, early_seq = bb.jd -fopn.jd)

ggplot(offset2, aes(year,buds_seq))+ geom_point()+(geom_abline(intercept=0, slope=0))+facet_wrap(~species)+ggtitle("Ind Leaf vs. Flower Budbust")+geom_smooth(method = "lm", se = FALSE)
ggplot(offset2, aes(year,functional_seq))+ geom_point()+(geom_abline(intercept=0, slope=0))+facet_wrap(~species)+ggtitle("Ind Leaf vs. Flower Developement")+geom_smooth(method = "lm", se = FALSE)
ggplot(offset2, aes(year,early_seq))+ geom_point()+(geom_abline(intercept=0, slope=0))+facet_wrap(~species)+ggtitle("Ind Leaf Budburst vs. Flower Developement")+geom_smooth(method = "lm", se = FALSE)
####It would be cool to over lay these plots with climate data

#######################################################################################################
#Now lets look at growing degree differences
######################################################################################################
##adding climate data
hf<- read.csv("data/WeatherData.csv", header=TRUE)
hf<-distinct(hf,Date,.keep_all = TRUE)

###calculating growing degree days (Cat's code with a little of my trouble shooting)
hf$gdd <- hf$AirT - 5
hf$gdd <-ifelse(hf$gdd>0, hf$gdd, 0)
hf$gdd <-ifelse(!is.na(hf$gdd), hf$gdd, 0) #added by dan-this is the problem line
hf$count <- ave(
  hf$gdd, hf$Year, 
  FUN=function(x) cumsum(c(0, head(x, -1)))
)
##clean to merge into one dataset
hf<-filter(hf,Year>=1990)
library(plyr)
hf<-rename(hf, c("Year"="year"))

d<- as.data.frame(rapply(object = d, f = round, classes = "numeric", how = "replace", digits = 0)) 
#creat combined dtata frame
df<-left_join(hf, d)

###plot Buds GDD
df2<-df%>%
  dplyr::select(year,species, JD, bb.jd, tree.id,count)
df2$day<- ifelse(df2$JD==df2$bb.jd,df2$JD,NA)
df2<-na.omit(df2)
ggplot(df2, aes(x=year, y=count)) + geom_point()+ggtitle("bud burst")+facet_wrap(~species)

####plot Flobud GDD
df3<-df%>%
  dplyr::select(year,species, JD, fbb.jd, tree.id,count)
df3$day<- ifelse(df3$JD==df3$fbb.jd,df3$JD,NA)
df3<-na.omit(df3)
ggplot(df3, aes(x=year, y=count)) + geom_point(aes(col="pink"))+ggtitle("floral bud burst")+facet_wrap(~species)
###extract standard deviation
##the function from my old srcipt is behaving badly. going to try again later
recap.bb<-df2 %>% dplyr::group_by(tree.id,species) %>% dplyr::summarise(mean.bb = mean(count),sd.bb=sd(count))
recap.fbb<-df3 %>% dplyr::group_by(tree.id,species) %>% dplyr::summarise(mean.fbb = mean(count),sd.fbb=sd(count))
##merge these sets
GDD.bud.com<-left_join(recap.bb, recap.fbb)
GDD.bud.com<-dplyr::select(GDD.bud.com,-mean.fbb)
GDD.bud.com<-gather(GDD.bud.com,phenophase,dev,sd.bb:sd.fbb)
ggplot(GDD.bud.com, aes(species, dev, color = phenophase))+geom_point()
GDD.bud.com<-filter(GDD.bud.com, species %in% c( "ACPE","ACRU", "ACSA","BEAL","FRAM","QURU"))
ggplot(GDD.bud.com, aes(species, dev, color=phenophase))+geom_point()
