####This is the master Rcode for Dan B's hysteranthy paper

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")

###Data#####################################################################
#data1
aln<-read.csv("datasheets_derived/alnus_delta_hyst.csv",header=TRUE)
frax<-read.csv("datasheets_derived/fraxinus_delta_hyst.csv",header=TRUE)
aes<-read.csv("datasheets_derived/aes_delta_hyst.csv",header=TRUE)

###analysis 1: Is hysteranthy changing with time? featureing PEP725 data

###function for ordering pepsite
pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)
  
}

###making the hinge
aln <- pepnumber(aln, "s_id")

aln$YEAR.hin <- aln$year
aln$YEAR.hin[which(aln$YEAR.hin<1980)] <- 1980
aln$YEAR.hin <- aln$YEAR.hin-1980

frax <- pepnumber(frax, "s_id")
frax$YEAR.hin <- frax$year
frax$YEAR.hin[which(frax$YEAR.hin<1980)] <- 1980
frax$YEAR.hin <- frax$YEAR.hin-1980

aes<- pepnumber(aes, "s_id")
aes$YEAR.hin <- aes$year
aes$YEAR.hin[which(aes$YEAR.hin<1980)] <- 1980
aes$YEAR.hin <- aes$YEAR.hin-1980

##### Alnus model
fit.aln.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aln) 
aln.dat<-as.data.frame(coef(fit.aln.brms))

aln.dat<-rownames_to_column(aln.dat, var = "peporder")
aln.dat<-merge(aln.dat,aln)
colnames(aln.dat)

names(aln.dat)[1]<-"peporder"
names(aln.dat)[4]<-"Intercept.2.5"
names(aln.dat)[5]<-"Intercept.97.5"
names(aln.dat)[2]<-"Intercept"
names(aln.dat)[6]<-"slope"
names(aln.dat)[7]<-"slope.2.5"
names(aln.dat)[8]<-"slope.97.5"
aln.dat<-dplyr::select(aln.dat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
aln.dat$hinge<-ifelse(aln.dat$YEAR.hin>0,1,0)

alphaALN<-mean(aln.dat$Intercept)
alphaALNlow<-mean(aln.dat$Intercept.2.5)
alphaALNhigh<-mean(aln.dat$Intercept.97.5)
betaALN<-mean(aln.dat$slope)
#betaALNlow<-mean(aln.dat$slope.2.5)
#betaALNhigh<-mean(aln.dat$slope.97.5)

####Fraxinus model
fit.frax.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=frax) 
frax.dat<-as.data.frame(coef(fit.frax.brms))
frax.dat<-rownames_to_column(frax.dat, var = "peporder")
frax.dat<-merge(frax.dat,frax)

names(frax.dat)[1]<-"peporder"
names(frax.dat)[4]<-"Intercept.2.5"
names(frax.dat)[5]<-"Intercept.97.5"
names(frax.dat)[2]<-"Intercept"
names(frax.dat)[6]<-"slope"
names(frax.dat)[7]<-"slope.2.5"
names(frax.dat)[8]<-"slope.97.5"
frax.dat<-dplyr::select(frax.dat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
frax.dat$hinge<-ifelse(frax.dat$YEAR.hin>0,1,0)


alphaFRAX<-mean(frax.dat$Intercept)
alphaFRAXlow<-mean(frax.dat$Intercept.2.5)
alphaFRAXhigh<-mean(frax.dat$Intercept.97.5)
betaFRAX<-mean(frax.dat$slope)


###Aesculus model
fit.aes.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aes) 
aes.dat<-as.data.frame(coef(fit.aes.brms))
aes.dat<-rownames_to_column(aes.dat, var = "peporder")
aes.dat<-merge(aes.dat,aes)

names(aes.dat)[1]<-"peporder"
names(aes.dat)[4]<-"Intercept.2.5"
names(aes.dat)[5]<-"Intercept.97.5"
names(aes.dat)[2]<-"Intercept"
names(aes.dat)[6]<-"slope"

names(aes.dat)[7]<-"slope.2.5"
names(aes.dat)[8]<-"slope.97.5"
aes.dat<-dplyr::select(aes.dat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
aes.dat$hinge<-ifelse(aes.dat$YEAR.hin>0,1,0)


alphaAES<-mean(aes.dat$Intercept)
betaAES<-mean(aes.dat$slope)
alphaAESlow<-mean(aes.dat$Intercept.2.5)
alphaAEShigh<-mean(aes.dat$Intercept.97.5)



allbetas<-aln.dat$slope
allalphas<-aln.dat$Intercept
allbetas1<-frax.dat$slope
allalphas1<-frax.dat$Intercept
allbetas2<-aes.dat$slope
allalphas2<-aes.dat$Intercept

# graph

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
segments(x0=1960,y0=alphaALN, x1=1980,y1=alphaALN,col="darkgreen",lwd=3)
segments(x0=1980,y0=alphaALN, x1=2015,y1=alphaALN+betaALN*35,col="darkgreen",lwd=3)
segments(x0=1960,y0=alphaALNlow, x1=2015,y1=alphaALNlow, lty=2, col="darkgreen",lwd=2)
segments(x0=1960,y0=alphaALNhigh, x1=2015,y1=alphaALNhigh, lty=2, col="darkgreen",lwd=2)
segments(x0=1960,y0=alphaFRAX, x1=1980,y1=alphaFRAX,col="red",lwd=3)
segments(x0=1980,y0=alphaFRAX, x1=2015,y1=alphaFRAX+betaFRAX*35,col="red",lwd=3)
segments(x0=1960,y0=alphaFRAXlow, x1=2015,y1=alphaFRAXlow, lty=2,col="red",lwd=2)
segments(x0=1960,y0=alphaFRAXhigh, x1=2015,y1=alphaFRAXhigh, lty=2,col="red",lwd=2)
segments(x0=1960,y0=alphaAES, x1=1980,y1=alphaAES,col="blue",lwd=3)
segments(x0=1980,y0=alphaAES, x1=2015,y1=alphaAES+betaAES*35,col="blue",lwd=3)
segments(x0=1960,y0=alphaAESlow, x1=2015,y1=alphaAESlow, lty=2,col="blue",lwd=2)
segments(x0=1960,y0=alphaAEShigh, x1=2015,y1=alphaAEShigh, lty=2,col="blue",lwd=2)
legend(2006,3, legend=c("A. glutinosa", "F. excelsior","A. hippocastanum"),col=c("darkgreen", "red","blue"), lwd=2, cex=0.6)
dev.off()


#####harvardforeest interannaulo

save.image("paper_full_analysis.RData")
