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
#jpeg("FLS_climate_change.jpeg",width = 6, height = 5, units = 'in', res = 300)
tiff("FLS_climate_change.tiff",width = 6, height = 5, units = 'in', res = 300)
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

###addd box plots for each species
all.data<-rbind(aesc50,alnu50,betu50,frax50)

prehinge<-filter(all.data,YEAR.hin==0)
recenty<-filter(all.data,YEAR.hin %in% c(30:35))
recenty$period<-"2010-2015"
prehinge$period<-"pre-1980"
hingy<-rbind(prehinge,recenty)
hingy$period <- factor(hingy$period, levels = c("pre-1980","2010-2015"))
scale_x_discrete("",labels=c("A. hippocastenum", "A. glutinosa","B. pendula","F. excelsior"))

setEPS()
postscript("pepboxplots.eps",width = 6, height = 3)
aeshingy<-filter(hingy,plant_id==101)

upper.limit <- quantile(aeshingy$FLS)[4] + 1.6*IQR(aeshingy$FLS)
lower.limit <- quantile(aeshingy$FLS)[2] - 1.6*IQR(aeshingy$FLS)
?geom_boxplot()

setEPS()
postscript("pepboxplots.aes.eps",width = 3, height = 5)
ggplot(aeshingy,aes(period,FLS))+geom_boxplot(aes(),outlier.shape = NA)+ggthemes::theme_base(base_size = 10)+
  scale_color_manual(values=c("grey68","black"))+geom_hline(yintercept=0, linetype="dotted")+
  scale_x_discrete("A. hippocastenum")+
 theme(axis.text.x = element_text(hjust = .5, vjust = .5),axis.title.x = element_text(face = "italic") )+coord_cartesian(ylim=c(lower.limit, upper.limit))
dev.off()


alhingy<-filter(hingy,plant_id==102)

upper.limit <- quantile(alhingy$FLS)[4] + 1.5*IQR(alhingy$FLS)
lower.limit <- quantile(alhingy$FLS)[2] - 1.5*IQR(alhingy$FLS)

setEPS()
postscript("pepboxplots.al.eps",width = 3, height = 5)
ggplot(alhingy,aes(period,FLS))+geom_boxplot(aes(),outlier.shape = NA)+ggthemes::theme_base(base_size = 10)+
  scale_color_manual(values=c("grey68","black"))+geom_hline(yintercept=0, linetype="dotted")+
  scale_x_discrete("A. glutinosa")+
  theme(axis.text.x = element_text(hjust = .5, vjust = .5),axis.title.x = element_text(face = "italic") )+coord_cartesian(ylim=c(lower.limit, upper.limit))
dev.off()


fraxhingy<-filter(hingy,plant_id==120)

upper.limit <- quantile(fraxhingy$FLS)[4] + 1.5*IQR(fraxhingy$FLS)
lower.limit <- quantile(fraxhingy$FLS)[2] - 1.5*IQR(fraxhingy$FLS)

setEPS()
postscript("pepboxplots.fr.eps",width = 3, height = 5)
ggplot(fraxhingy,aes(period,FLS))+geom_boxplot(aes(),outlier.shape = NA)+ggthemes::theme_base(base_size = 10)+
  scale_color_manual(values=c("grey68","black"))+geom_hline(yintercept=0, linetype="dotted")+
  scale_x_discrete("F. excelsior")+
  theme(axis.text.x = element_text(hjust = .5, vjust = .5),axis.title.x = element_text(face = "italic") )+coord_cartesian(ylim=c(lower.limit, upper.limit))
dev.off()


bethingy<-filter(hingy,plant_id==106)

upper.limit <- quantile(bethingy$FLS)[4] + 1.5*IQR(bethingy$FLS)
lower.limit <- quantile(bethingy$FLS)[2] - 1.5*IQR(bethingy$FLS)

setEPS()
postscript("pepboxplots.bet.eps",width = 3, height = 5)
ggplot(bethingy,aes(period,FLS))+geom_boxplot(aes(),outlier.shape = NA)+ggthemes::theme_base(base_size = 10)+
  scale_color_manual(values=c("grey68","black"))+geom_hline(yintercept=0, linetype="dotted")+
  scale_x_discrete("B. pendula")+
  theme(axis.text.x = element_text(hjust = .5, vjust = .5),axis.title.x = element_text(face = "italic") )+coord_cartesian(ylim=c(lower.limit, upper.limit))
dev.off()



setEPS()
postscript("FLS_climate_change.eps",width = 8, height = 7)
plot(c(1960,2015), c(-20,30), type = "n", xlab = "Year", ylab = "Days between flowering and leafing", bty='l')
rect(xleft=1960,xright=2015,ybottom=fixef(fit.aln50.brms)[1,3],ytop=fixef(fit.aln50.brms)[1,4],col="tomato",border=NA)
segments(x0=1960,y0=fixef(fit.aln50.brms)[1], x1=1980, y1=fixef(fit.aln50.brms)[1],col="red3",lwd=3)
segments(x0=1980,y0=fixef(fit.aln50.brms)[1], x1=2015,y1=fixef(fit.aln50.brms)[1]+fixef(fit.aln50.brms)[2]*35,col="red3",lwd=3)


rect(xleft=1960,xright=2015,ybottom=fixef(fit.frax50.brms)[1,3],ytop=fixef(fit.frax50.brms)[1,4],col="khaki1",border=NA)
segments(x0=1960,y0=fixef(fit.frax50.brms)[1], x1=1980,y1=fixef(fit.frax50.brms)[1],col="darkgoldenrod1",lwd=3)
segments(x0=1980,y0=fixef(fit.frax50.brms)[1], x1=2015,y1=fixef(fit.frax50.brms)[1]+fixef(fit.frax50.brms)[2]*35,col="darkgoldenrod1",lwd=3)

rect(xleft=1960,xright=2015,ybottom=fixef(fit.aesc50.brms)[1,3],ytop=fixef(fit.aesc50.brms)[1,4],col="lightblue",border=NA)
segments(x0=1960,y0=fixef(fit.aesc50.brms)[1], x1=1980,y1=fixef(fit.aesc50.brms)[1],col="navyblue",lwd=3)
segments(x0=1980,y0=fixef(fit.aesc50.brms)[1], x1=2015,y1=fixef(fit.aesc50.brms)[1]+fixef(fit.aesc50.brms)[2]*35,col="navyblue",lwd=3)

rect(xleft=1960,xright=2015,ybottom=fixef(fit.betu50.brms)[1,3],ytop=fixef(fit.betu50.brms)[1,4],col="springgreen3",border=NA)
segments(x0=1960,y0=fixef(fit.betu50.brms)[1], x1=1980,y1=fixef(fit.betu50.brms)[1],col="darkgreen",lwd=3)
segments(x0=1980,y0=fixef(fit.betu50.brms)[1], x1=2015,y1=fixef(fit.betu50.brms)[1]+fixef(fit.betu50.brms)[2]*35,col="darkgreen",lwd=3)

text(x=1965, 24,font=3,cex=.8,"Alnus glutinosa")
text(x=1966, 12,font=3,cex=.8,"Fraxinus excelsior")
text(x=1965, 4,font=3,cex=.8,"Betula pendula")
text(x=1969, -16,font=3,cex=.8,"Aesculus hippocastanum")
dev.off()
