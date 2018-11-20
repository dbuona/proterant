####brms hinge models

####only fraxinus, alnus, fagus, tilia and aesculus data has been updated using agood standard (stations with 50 or more years of data)
####if you want to run other models, fix the data on delta_hyst.R first. But only alnus, ausculus and frax have more than 1 or 2 stations
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
load("hinges.Rdata")
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

#read in clean data
aln<-read.csv("alnus_delta_hyst.csv",header=TRUE)
frax<-read.csv("fraxinus_delta_hyst.csv",header=TRUE)
aes<-read.csv("aes_delta_hyst.csv",header=TRUE)
#fag<-read.csv("fag_delta_hyst.csv",header=TRUE)
#til<-read.csv("tilia_delta_hyst.csv",header=TRUE)

###function for ordering pepsite
pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)
  
}

###making the hinge
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

##### running the models
fit.aln.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aln) 
dat<-as.data.frame(coef(fit.aln.brms))

dat<-rownames_to_column(dat, var = "peporder")
newdat<-merge(dat,aln)
colnames(newdat)

names(newdat)[1]<-"peporder"
names(newdat)[4]<-"Intercept.2.5"
names(newdat)[5]<-"Intercept.97.5"
names(newdat)[2]<-"Intercept"
names(newdat)[6]<-"slope"
names(newdat)[7]<-"slope.2.5"
names(newdat)[8]<-"slope.97.5"
newdata<-dplyr::select(newdat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
newdata$hinge<-ifelse(newdata$YEAR.hin>0,1,0)
colnames(newdata)
alphaX<-mean(newdata$Intercept)
alphaXlow<-mean(newdata$Intercept.2.5)
alphaXhigh<-mean(newdata$Intercept.97.5)
betaX<-mean(newdata$slope)
betaXlow<-mean(newdata$slope.2.5)
betaXhigh<-mean(newdata$slope.97.5)

allbetas<-newdata$slope
allalphas<-newdata$Intercept
allbetas1<-newdata1$slope
allalphas1<-newdata1$Intercept
allbetas2<-newdata2$slope
allalphas2<-newdata2$Intercept

# graph
dev.off()
jpeg("..//figure/FLS_climate_change.jpeg")
plot(c(1960,2015), c(-30,45), type = "n", xlab = "year", ylab = "FLS offset", bty='l')
#rect(xleft=1960, ybottom=-30, xright=1980, ytop=50,col="ivory1" )
#rect(xleft=1980, ybottom=-30, xright=2015, ytop=50,col="ivory2")
segments(x0=1980,y0=allalphas,x1=2015,y1=allalphas+allbetas*35,col="azure3",lty="dotted" )
segments(x0=1960,y0=allalphas,x1=1980,y1=allalphas,col="azure3",lty="dotted")
segments(x0=1980,y0=allalphas1,x1=2015,y1=allalphas1+allbetas1*35,col="azure3")
segments(x0=1960,y0=allalphas1,x1=1980,y1=allalphas1,col="azure3")
segments(x0=1980,y0=allalphas2,x1=2015,y1=allalphas2+allbetas2*35,col="azure3",lty="dashed")
segments(x0=1960,y0=allalphas2,x1=1980,y1=allalphas2,col="azure3",lty="dashed")
segments(x0=1960,y0=alphaX, x1=1980,y1=alphaX,col="darkgreen",lwd=3)
segments(x0=1980,y0=alphaX, x1=2015,y1=alphaX+betaX*35,col="darkgreen",lwd=3)
segments(x0=1960,y0=alphaXlow, x1=2015,y1=alphaXlow, lty=2, col="darkgreen",lwd=2)
segments(x0=1960,y0=alphaXhigh, x1=2015,y1=alphaXhigh, lty=2, col="darkgreen",lwd=2)
segments(x0=1960,y0=alpha1, x1=1980,y1=alpha1,col="red",lwd=3)
segments(x0=1980,y0=alpha1, x1=2015,y1=alpha1+beta1*35,col="red",lwd=3)
segments(x0=1960,y0=alphalow, x1=2015,y1=alphalow, lty=2,col="red",lwd=2)
segments(x0=1960,y0=alphahigh, x1=2015,y1=alphahigh, lty=2,col="red",lwd=2)
segments(x0=1960,y0=alpha2, x1=1980,y1=alpha2,col="blue",lwd=3)
segments(x0=1980,y0=alpha2, x1=2015,y1=alpha2+beta2*35,col="blue",lwd=3)
segments(x0=1960,y0=alpha2low, x1=2015,y1=alpha2low, lty=2,col="blue",lwd=2)
segments(x0=1960,y0=alpha2high, x1=2015,y1=alpha2high, lty=2,col="blue",lwd=2)
legend(2002,3, legend=c("A. glutinosa", "F. excelsior","A. hippocastanum"),col=c("darkgreen", "red","blue"), lwd=2, cex=0.6)

dev.off()




plottyX<-ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980),size=0.1,color="lightskyblue3",linetype="solid")+geom_segment(aes(y=alphaX,yend=alphaX,x=1960, xend=1980),color="slategray",linetype="solid",size=1)+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope),size=0.1,color="lightsalmon1",linetype="solid",linejoin = c('bevel'))+geom_segment(aes(x=1980,xend=2015,y=alphaX,yend=alphaX+(35*betaX)),color="orchid4",linetype="solid",size=1)+theme_linedraw()+theme(legend.position="none")+ggtitle("Alnus glutinosa")
plottyalnus<-plottyX+geom_segment(aes(y=alphaXlow,yend=alphaXlow,x=1960, xend=1980),data=newdat,color="blue",linetype="dashed",size=0.5)+geom_segment(aes(y=alphaXhigh,yend=alphaXhigh, x=1960, xend=1980),data=newdat,color="blue",linetype="dashed",size=0.5)

alnusplot.close<-plottyalnus+geom_segment(aes(y=alphaXlow,yend=alphaXlow,x=1980, xend=2015),data=newdat,color="black",linetype="dashed",size=0.5)+geom_segment(aes(y=alphaXhigh,yend=alphaXhigh, x=1980, xend=2015),data=newdat,color="black",linetype="dashed",size=0.5)
alnus.plotfinal<-alnusplot.close+geom_segment(aes(x=1980,xend=2015,y=alphaXhigh,yend=alphaXhigh+35*betaXhigh))

####Fraxinus model
fit.frax.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=frax) 
dat1<-as.data.frame(coef(fit.frax.brms))
dat1<-rownames_to_column(dat1, var = "peporder")
newdat1<-merge(dat1,frax)
colnames(newdat1)
names(newdat1)[1]<-"peporder"
names(newdat1)[4]<-"Intercept.2.5"
names(newdat1)[5]<-"Intercept.97.5"
names(newdat1)[2]<-"Intercept"
names(newdat1)[6]<-"slope"
names(newdat1)[7]<-"slope.2.5"
names(newdat1)[8]<-"slope.97.5"
newdata1<-dplyr::select(newdat1,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
newdata1$hinge<-ifelse(newdata1$YEAR.hin>0,1,0)


alpha1<-mean(newdata1$Intercept)
alphalow<-mean(newdata1$Intercept.2.5)
alphahigh<-mean(newdata1$Intercept.97.5)
beta1<-mean(newdata1$slope)


#ggplot(newdata1,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha1,yend=alpha1,x=1960, xend=1980),color="red")
#ggplot(newdata1,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha1,yend=alpha1+35*beta1),color="red")

plotty2<-ggplot(newdata1,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980),size=0.1,color="lightgray")+geom_segment(aes(y=alpha1,yend=alpha1,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope),size=0.1,color="lightgray")+geom_segment(aes(x=1980,xend=2015,y=alpha1,yend=alpha1+35*beta1),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Fraxinus excelsior")
plottyfrax<-plotty2+geom_segment(aes(y=alphalow,yend=alphalow,x=1960, xend=2015),data=newdat1,color="blue",linetype="dashed",size=0.5)+geom_segment(aes(y=alphahigh,yend=alphahigh,x=1960, xend=2015),data=newdat1,color="blue",linetype="dashed",size=0.5)

###Aesculus model

fit.aes.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aes) 
dat2<-as.data.frame(coef(fit.aes.brms))
dat2<-rownames_to_column(dat2, var = "peporder")
newdat2<-merge(dat2,aes)
colnames(newdat2)
names(newdat2)[1]<-"peporder"
names(newdat2)[4]<-"Intercept.2.5"
names(newdat2)[5]<-"Intercept.97.5"
names(newdat2)[2]<-"Intercept"
names(newdat2)[6]<-"slope"

names(newdat2)[7]<-"slope.2.5"
names(newdat2)[8]<-"slope.97.5"
newdata2<-dplyr::select(newdat2,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)


newdata2$hinge<-ifelse(newdata2$YEAR.hin>0,1,0)


alpha2<-mean(newdata2$Intercept)
beta2<-mean(newdata2$slope)
alpha2low<-mean(newdata2$Intercept.2.5)
alpha2high<-mean(newdata2$Intercept.97.5)


plotty3<-ggplot(newdata2,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980),size=0.1,color="lightgray")+geom_segment(aes(y=alpha2,yend=alpha2,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope),size=0.1,color="lightgray")+geom_segment(aes(x=1980,xend=2015,y=alpha2,yend=alpha2+35*beta2),color="red")+theme_tufte()+theme(legend.position="none")+ggtitle("Aesculus hippocastanum")
plottyaes<-plotty3+geom_segment(aes(y=alpha2low,yend=alpha2low,x=1960, xend=2015),data=newdat2,color="blue",linetype="dashed",size=0.5)+geom_segment(aes(y=alpha2high,yend=alpha2high,x=1960, xend=2015),data=newdat2,color="blue",linetype="dashed",size=0.5)

?geom_segment()

save.image(file="hinges.Rdata")

### are there any clines
FRAXINUS.DATA<-newdata1
unique(FRAXINUS.DATA$lat)
quantile(FRAXINUS.DATA$lat)
median(FRAXINUS.DATA$lat)



###STOP HERE FOR NOW--maps below arent usable and other species don't have much data#########################
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




library(gridExtra)
Alna<-grid.arrange(plotty,mappy,ncol=2, nrow=1)
fraxer<-grid.arrange(plottyfrax,mappy1,ncol=2, nrow=1)
aescu<-grid.arrange(plottyaes,mappy2,ncol=2, nrow=1)
fagus<-grid.arrange(plotty4,mappy3,ncol=2, nrow=1)
tilia<-grid.arrange(plotty5,mappy4,ncol=2, nrow=1)

save.image(file="hinges.Rdata")
