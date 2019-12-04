rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

setwd("~/Documents/git/proterant/sub_projs/PEP725/")
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(tibble)
load("..//PEP725/hinge.models")

frax50<-read.csv("input/frax50year.csv")
aesc50<-read.csv("input/aesc50year.csv")
alnu50<-read.csv("input/alnu50year.csv")
betu50<-read.csv("input/betu50year.csv")

###function for ordering pepsite
pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)
  
}

makehinge<-function(x){
  x<- pepnumber(x, "s_id") 
x<-transform(x,YEAR.hin=ifelse(year<=1980,1980,year))
x<-transform(x,YEAR.hin=YEAR.hin-1980)}

alnu50<-makehinge(alnu50)
aesc50<-makehinge(aesc50)
betu50<-makehinge(betu50)
frax50<-makehinge(frax50)

##now model each species separately
fit.aln50.brms<-brm(FLS~YEAR.hin+(YEAR.hin|peporder),data=alnu50) 
fit.betu50.brms<-brm(FLS~YEAR.hin+(YEAR.hin|peporder),data=betu50) 
fit.frax50.brms<-brm(FLS~YEAR.hin+(YEAR.hin|peporder),data=frax50) 
fit.aesc50.brms<-brm(FLS~YEAR.hin+(YEAR.hin|peporder),data=aesc50)

fixef(fit.frax50.brms)[1,2]
jpeg("FLS_climate_change.jpeg",width = 6, height = 5, units = 'in', res = 300)
plot(c(1960,2015), c(-20,30), type = "n", xlab = "Year", ylab = "Days between flowering and leafing", bty='l')
segments(x0=1960,y0=fixef(fit.aln50.brms)[1], x1=1980, y1=fixef(fit.aln50.brms)[1],col="red",lwd=3)
segments(x0=1980,y0=fixef(fit.aln50.brms)[1], x1=2015,y1=fixef(fit.aln50.brms)[1]+fixef(fit.aln50.brms)[2]*35,col="red",lwd=3)
segments(x0=1960,y0=coef(fit.aln50.brms)[1], x1=1980, y1=coef(fit.aln50.brms)[1],col="red",lwd=3)
rect(xleft=1960,xright=2015,ybottom=fixef(fit.aln50.brms)[1,3],ytop=fixef(fit.aln50.brms)[1,4],col=rgb(1,0,0,alpha=0.4),border=NA)

segments(x0=1960,y0=fixef(fit.frax50.brms)[1], x1=1980,y1=fixef(fit.frax50.brms)[1],col="darkgoldenrod1",lwd=3)
segments(x0=1980,y0=fixef(fit.frax50.brms)[1], x1=2015,y1=fixef(fit.frax50.brms)[1]+fixef(fit.frax50.brms)[2]*35,col="darkgoldenrod1",lwd=3)
rect(xleft=1960,xright=2015,ybottom=fixef(fit.frax50.brms)[1,3],ytop=fixef(fit.frax50.brms)[1,4],col=rgb(1,.9,0,alpha=0.4),border=NA)

segments(x0=1960,y0=fixef(fit.aesc50.brms)[1], x1=1980,y1=fixef(fit.aesc50.brms)[1],col="navyblue",lwd=3)
segments(x0=1980,y0=fixef(fit.aesc50.brms)[1], x1=2015,y1=fixef(fit.aesc50.brms)[1]+fixef(fit.aesc50.brms)[2]*35,col="navyblue",lwd=3)
rect(xleft=1960,xright=2015,ybottom=fixef(fit.aesc50.brms)[1,3],ytop=fixef(fit.aesc50.brms)[1,4],col=rgb(.4,.4,1,alpha=0.3),border=NA)

segments(x0=1960,y0=fixef(fit.betu50.brms)[1], x1=1980,y1=fixef(fit.betu50.brms)[1],col="darkgreen",lwd=3)
segments(x0=1980,y0=fixef(fit.betu50.brms)[1], x1=2015,y1=fixef(fit.betu50.brms)[1]+fixef(fit.betu50.brms)[2]*35,col="darkgreen",lwd=3)
rect(xleft=1960,xright=2015,ybottom=fixef(fit.betu50.brms)[1,3],ytop=fixef(fit.betu50.brms)[1,4],col=rgb(.3,1,.1,alpha=0.3),border=NA)

text(x=1965, 24,font=3,cex=.8,"Alnus glutinosa")
text(x=1966, 12,font=3,cex=.8,"Fraxinus excelsior")
text(x=1965, 4,font=3,cex=.8,"Betula pendula")
text(x=1969, -16,font=3,cex=.8,"Aesculus hippocastanum")
dev.off()
save.image("..//PEP725/hinge.models")


