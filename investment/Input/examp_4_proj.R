## 1) house keeping

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

### working direct
setwd("~/Documents/git/proterant/investment/Input")

## 2 packages
#install.packages("dplyr")
library(dplyr)
library(ggplot2)

d<-read.csv("pruno_checker - pruno_checker.csv")


head(d) ## see the first road
nrow(d)

### Q1 how many flower did we have per dataset
table(d$flowers)
75/nrow(d)





#what percentage of specimens have flowers or fruits
length(which(d$flowers=="Y"))/nrow(d)

length(which(d$bbch.f %in% c(60:67)))/nrow(d)


sumdat<-as.data.frame(table(d$specificEpithet,d$flowers)[,2]/10)



d2<-na.omit(d)

ggplot(d2, aes(x = as.factor(bbch.v),as.factor( bbch.f)))+
   stat_bin2d(bins=c(7,6))+
geom_rect(xmin=as.factor(0),xmax=as.factor(17),ymin=as.factor(60),ymax=as.factor(65),alpha=0.1,fill="goldenrod")+
   facet_wrap(~specificEpithet)

+
   theme_minimal()
   
?stat_bin2d()

 d2$FLS<-ifelse(d2$bbch.f<67 &d2$bbch.v<=17,1,0)
ggplot(d2,aes(as.factor(FLS)))+geom_bar()+facet_wrap(~specificEpithet)

MOD1<-brms::brm(FLS~1+(1|specificEpithet),family="bernoulli",data=d2)

newdata<-data.frame(specificEpithet=unique(d2$specificEpithet))
pred<-fitted(MOD1,newdata=newdata)
newdata<-cbind(newdata,pred)
ggplot(newdata,aes(specificEpithet,Estimate))+geom_point()+geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0.1)
       