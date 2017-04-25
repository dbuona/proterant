##A sbuset of flo_v_leaf_master to make a figure for CSEE2017 on 4/21/17

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
d.1<-filter(d.spp, species %in% c( "ACPE","ACRU", "ACSA","BEAL","ILVE","NEMU","POTR","FRAM","QURU")) #filter species where observation go full legnth of study
d.2<-gather(d.1,phenophase, eventday, bb.jd:fopn.jd)

### rename everything
d.2$species[d.2$species == "ACPE"] <- "A. pennsylvanicum"
d.2$species[d.2$species == "ACRU"] <- "A. rubrum"
d.2$species[d.2$species == "ACSA"] <- "A. saccharrum"
d.2$species[d.2$species == "BEAL"] <- "B. allegheniensis"
d.2$species[d.2$species == "ILVE"] <- "I. verticillata"
d.2$species[d.2$species == "NEMU"] <- "I. mucronata"
d.2$species[d.2$species == "POTR"] <- "P. tremuloides"
d.2$species[d.2$species == "FRAM"] <- "F. americana"
d.2$species[d.2$species == "QURU"] <- "Q. rubra"
d.2$species[d.2$species == "ACPE"] <- "A.pennsylvanicum"
d.2$phenophase[d.2$phenophase == "bb.jd"] <- "budburst"
d.2$phenophase[d.2$phenophase == "fbb.jd"] <- "floral budburst"
d.2$phenophase[d.2$phenophase == "l75.jd"] <- "75% leaf size"
d.2$phenophase[d.2$phenophase == "fopn.jd"] <- "open flower"

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
View(slopes.ind)
write.csv(slopes.ind, file = "HF_slopes.csv",row.names=FALSE)

###Now shrubs
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/git/buds/analyses/Flo_buds/input")
fb<-read.csv("..//input/Budburst By Day.csv", header = TRUE)
#make subset for the good shurbs leftover from Mag 7
goodshrubs<-filter(fb, sp %in% c( "CORCOR","ILEMUC", "PRUPEN"))

goodshrubs<-group_by(goodshrubs,sp,treatcode)
goodshrubs<-gather(goodshrubs,phenophase,eventday,fday:lday)
goodshrubs<-filter(goodshrubs, treatcode %in% c( "CL0","CS0", "WL0","WS0"))


goodshrubs$sp[goodshrubs$sp == "ILEMUC"] <- "I. mucronata"
goodshrubs$sp[goodshrubs$sp == "PRUPEN"] <- "P. pensylvanica"
goodshrubs$sp[goodshrubs$sp == "CORCOR"] <- "C. cornuta"

goodshrubs$treatcode[goodshrubs$treatcode == "WS0"] <- "WS"
goodshrubs$treatcode[goodshrubs$treatcode == "CS0"] <- "CS"
goodshrubs$treatcode[goodshrubs$treatcode == "WL0"] <- "WL"
goodshrubs$treatcode[goodshrubs$treatcode == "CL0"] <- "CL"

goodshrubs$phenophase[goodshrubs$phenophase == "fday"] <- "flowering"
goodshrubs$phenophase[goodshrubs$phenophase == "lday"] <- "leafing"

goodshrubs$warm[goodshrubs$warm == 20] <- "warm"
goodshrubs$warm[goodshrubs$warm == 15] <- "cool"

goodshrubs$photo[goodshrubs$photo == 12] <- "long"
goodshrubs$photo[goodshrubs$photo == 8] <- "short"

colnames(goodshrubs)[7] <- "temperature"
colnames(goodshrubs)[8] <- "photoperiod"


q<-ggplot(goodshrubs, aes(x=temperature, y=eventday, color=photoperiod, shape=phenophase)) +stat_summary()
 q+facet_wrap(~sp)+labs(x="temperature", y="days since initiation")
 
 q<-ggplot(goodshrubs, aes(x=treatcode, y=eventday, color=phenophase)) +stat_summary()
 q+facet_wrap(~sp)+labs(x="temperature", y="days since initiation")



