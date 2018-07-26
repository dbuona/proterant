####brms hinge models
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
#load("hinges.Rdata")
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)
library(arm)
library(shinystan)
library(ggthemes)
library(brms)
library(tibble)

aln<-read.csv("alnus_delta_hyst.csv",header=TRUE)
frax<-read.csv("fraxinus_delta_hyst.csv",header=TRUE)
aes<-read.csv("aes_delta_hyst.csv",header=TRUE)
fag<-read.csv("fag_delta_hyst.csv",header=TRUE)
til<-read.csv("tilia_delta_hyst.csv",header=TRUE)

pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)
  
}


alnUSW <- pepnumber(aln, "s_id")
aln<-alnUSW
aln$YEAR.hin <- aln$year
aln$YEAR.hin[which(aln$YEAR.hin<1980)] <- 1980
aln$YEAR.hin <- aln$YEAR.hin-1980

fraxUSW <- pepnumber(frax, "s_id")
frax<-fraxUSW
frax$YEAR.hin <- frax$year
frax$YEAR.hin[which(frax$YEAR.hin<1980)] <- 1980
frax$YEAR.hin <- frax$YEAR.hin-1980

aesUSW <- pepnumber(aes, "s_id")
aes<-aesUSW
aes$YEAR.hin <- aes$year
aes$YEAR.hin[which(aes$YEAR.hin<1980)] <- 1980
aes$YEAR.hin <- aes$YEAR.hin-1980

fagUSW <- pepnumber(fag, "s_id")
fag<-fagUSW
fag$YEAR.hin <- fag$year
fag$YEAR.hin[which(fag$YEAR.hin<1980)] <- 1980
fag$YEAR.hin <- fag$YEAR.hin-1980

tilUSW <- pepnumber(til, "s_id")
til<-tilUSW
til$YEAR.hin <- til$year
til$YEAR.hin[which(til$YEAR.hin<1980)] <- 1980
til$YEAR.hin <- til$YEAR.hin-1980


fit.aln.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aln) 
dat<-as.data.frame(coef(fit.aln.brms))

dat<-rownames_to_column(dat, var = "peporder")
newdat<-merge(dat,aln)
colnames(newdat)
names(newdat)[1]<-"peporder"
names(newdat)[2]<-"Intercept"
names(newdat)[6]<-"slope"
newdata<-dplyr::select(newdat,peporder,Intercept,slope,s_id,year,offset,YEAR.hin,lat,lon)
newdata$hinge<-ifelse(newdata$YEAR.hin>0,1,0)

alpha<-mean(newdata$Intercept)
beta<-mean(newdata$slope)

#ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha,yend=alpha,x=1960, xend=1980),color="red")
#ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha,yend=alpha+35*beta),color="red")

plotty<-ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha,yend=alpha,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha,yend=alpha+35*beta),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Alnus")

####put them on a map

####Fraxinus model
fit.frax.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=frax) 
dat1<-as.data.frame(coef(fit.frax.brms))
dat1<-rownames_to_column(dat1, var = "peporder")
newdat1<-merge(dat1,frax)
colnames(newdat1)
names(newdat1)[1]<-"peporder"
names(newdat1)[2]<-"Intercept"
names(newdat1)[6]<-"slope"
newdata1<-dplyr::select(newdat1,peporder,Intercept,slope,s_id,year,offset,YEAR.hin,lat,lon)
newdata1$hinge<-ifelse(newdata1$YEAR.hin>0,1,0)


alpha1<-mean(newdata1$Intercept)
beta1<-mean(newdata1$slope)

#ggplot(newdata1,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha1,yend=alpha1,x=1960, xend=1980),color="red")
#ggplot(newdata1,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha1,yend=alpha1+35*beta1),color="red")

plotty2<-ggplot(newdata1,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha1,yend=alpha1,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha1,yend=alpha1+35*beta1),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Fraxinus")

###Aesculus model

fit.aes.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aes) 
dat2<-as.data.frame(coef(fit.aes.brms))
dat2<-rownames_to_column(dat2, var = "peporder")
newdat2<-merge(dat2,aes)
colnames(newdat2)
names(newdat2)[1]<-"peporder"
names(newdat2)[2]<-"Intercept"
names(newdat2)[6]<-"slope"
newdata2<-dplyr::select(newdat2,peporder,Intercept,slope,s_id,year,offset,YEAR.hin,lat,lon)
newdata2$hinge<-ifelse(newdata2$YEAR.hin>0,1,0)


alpha2<-mean(newdata2$Intercept)
beta2<-mean(newdata2$slope)

#ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha2,yend=alpha2,x=1960, xend=1980),color="red")
ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha2,yend=alpha2+35*beta2),color="red")

plotty3<-ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha2,yend=alpha2,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha2,yend=alpha2+35*beta2),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Aesculus")

####Fagus model
fit.fag.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=fag) 
dat3<-as.data.frame(coef(fit.fag.brms))
dat3<-rownames_to_column(dat3, var = "peporder")
newdat3<-merge(dat3,fag)
colnames(newdat3)
names(newdat3)[1]<-"peporder"
names(newdat3)[2]<-"Intercept"
names(newdat3)[6]<-"slope"
newdata3<-dplyr::select(newdat3,peporder,Intercept,slope,s_id,year,offset,YEAR.hin,lat,lon)
newdata3$hinge<-ifelse(newdata3$YEAR.hin>0,1,0)


alpha3<-mean(newdata3$Intercept)
beta3<-mean(newdata3$slope)

#ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha2,yend=alpha2,x=1960, xend=1980),color="red")
#ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha2,yend=alpha2+35*beta2),color="red")
plotty4<-ggplot(newdata3,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha3,yend=alpha3,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha3,yend=alpha3+35*beta3),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Fagus")

#######Tilia
fit.til.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=til) 
dat4<-as.data.frame(coef(fit.til.brms))
dat4<-rownames_to_column(dat4, var = "peporder")
newdat4<-merge(dat4,til)
colnames(newdat4)
names(newdat4)[1]<-"peporder"
names(newdat4)[2]<-"Intercept"
names(newdat4)[6]<-"slope"
newdata4<-dplyr::select(newdat4,peporder,Intercept,slope,s_id,year,offset,YEAR.hin,lat,lon)
newdata4$hinge<-ifelse(newdata4$YEAR.hin>0,1,0)


alpha4<-mean(newdata4$Intercept)
beta4<-mean(newdata4$slope)

#ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha2,yend=alpha2,x=1960, xend=1980),color="red")
#ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha2,yend=alpha2+35*beta2),color="red")
plotty5<-ggplot(newdata4,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha4,yend=alpha4,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha4,yend=alpha4+35*beta4),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Tilia")


library("ggmap")
library(mapproj)
myLocation<-c(lon=9.8050,lat=51.9)
myMap<-get_map(location=myLocation,source="stamen",maptype="toner",zoom=5)
mappy<-ggmap(myMap)+geom_point(aes(x = lon, y = lat, color=slope), data = newdata ,alpha = .5, size = 1)+scale_colour_gradientn(colours = rainbow(8),limits = c(-1.5, 2))

mappy1<-ggmap(myMap)+geom_point(aes(x = lon, y = lat, color=slope), data = newdata1 ,alpha = .5, size = 1)+scale_colour_gradientn(colours = rainbow(8),limits = c(-1.5, 2))

mappy2<-ggmap(myMap)+geom_point(aes(x = lon, y = lat, color=slope), data = newdata2 ,alpha = .5, size = 1)+scale_colour_gradientn(colours = rainbow(8),limits = c(-1.5, 2))

mappy3<-ggmap(myMap)+geom_point(aes(x = lon, y = lat, color=slope), data = newdata3 ,alpha = .5, size = 1)+scale_colour_gradientn(colours = rainbow(8),limits = c(-1.5, 2))

mappy4<-ggmap(myMap)+geom_point(aes(x = lon, y = lat, color=slope), data = newdata4 ,alpha = .5, size = 1)+scale_colour_gradientn(colours = rainbow(8),limits = c(-1.5, 2))


library(gridExtra)
Alna<-grid.arrange(plotty,mappy,ncol=2, nrow=1)
fraxer<-grid.arrange(plotty2,mappy1,ncol=2, nrow=1)
aescu<-grid.arrange(plotty3,mappy2,ncol=2, nrow=1)
fagus<-grid.arrange(plotty4,mappy3,ncol=2, nrow=1)
tilia<-grid.arrange(plotty5,mappy4,ncol=2, nrow=1)

save.image(file="hinges.Rdata")
