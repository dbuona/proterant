### Explore the prunus data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(ggplot2)
library(brms)



setwd("~/Documents/git/proterant/investment/input")

flowers<-read.csv("flowermeasures - measurements.csv")
book<-read.csv("..//Data/zarchieve/cherry_data.csv")


pdsi<-read.csv("raw_means_pdsi.csv")
#fls<-read.csv("fls_gaus_estimate.csv")

meas<-flowers %>% group_by(specificEpithet) %>% dplyr::summarise(petaminl=min(pental_lengh_mm,na.rm=TRUE),pedtalmax=max(pental_lengh_mm,na.rm=TRUE),
                                                                 petalmean=mean(pental_lengh_mm,na.rm=TRUE),petalsd=sd(pental_lengh_mm,na.rm=TRUE),petalmed=median(pental_lengh_mm,na.rm=TRUE))
colnames(meas)[1]<-"species"

book<-filter(book,species %in% unique(meas$species))
book<-dplyr::select(book, 1:3)
book<-left_join(book,meas)
book$median<-(book$petal_low+book$petal_high)/2

book$species <- factor(book$species, levels = book$species[order(book$median)])
pd=position_dodge(width=0.8)

jpeg("../Plots/petal_length_comps.jpeg")
ggplot(book,aes(species,petal_low))+
  geom_linerange(aes(ymin=petaminl,ymax=pedtalmax),color="red",alpha=0.2,size=6)+
  geom_linerange(aes(ymin=petal_low,ymax=petal_high),size=4,alpha=0.4)+
  geom_point(aes(species,petalmean),size=4,color="red")+
  geom_linerange(aes(ymin=petalmean-petalsd,ymax=petalmean+petalsd),color="red")+ggthemes::theme_clean()+ylab("petal length")
dev.off()
colnames(pdsi)[1]<-"species"
pdsi<-filter(pdsi, species %in% unique(fls$species))


fls<-left_join(fls,pdsi)
fls<-left_join(fls,orde)

ggplot(fls,aes(Estimate,meanpdsi))+geom_point()+stat_smooth(method="lm")
ggplot(fls,aes(Estimate,minpdsi))+geom_point()+stat_smooth(method="lm")
ggplot(fls,aes(Estimate,petal))+geom_point()+stat_smooth(method="lm")

ggplot(fls,aes(meanpdsi,petal))+geom_point()+stat_smooth(method="lm")
ggplot(fls,aes(minpdsi,petal))+geom_point()+stat_smooth(method="lm")
ggplot(fls,aes(minpdsi,meanpdsi))+geom_point()+stat_smooth(method="lm")

get_prior(petal|mi(pedtalsd)~me(Estimate,SE),data=fls)


mod1<-brm(minpdsi~me(Estimate,SE),data=fls)

#zscore everything
function(x){}
mean(x)
lonpoints<-d$lon # make vector of prunus coordinates
latpoints<-d$lat #