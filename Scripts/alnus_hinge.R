rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")

load(file="stan/fit.hinge.aln.Rda") #skip to line 103
####cant get stan code to run
runstan = TRUE
mapsandsummaries = TRUE # these are just summaries for mapping (they're slow but not terribly), note that right now they run for the 10+ year dataset
use20yrs = FALSE # Otherwise it uses all data with 10 or more years
shinystancheck = FALSE # If you want to look at output in shiny stan
ncores = 2 # how many cores to use, only applies when runstan=TRUE

## Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(rgdal)
library(ggplot2)
library(lubridate)
library(rstan)
library(arm)
library(shinystan)
library(ggthemes)
aln<-read.csv("alnus_delta_hyst.csv",header=TRUE)


alnALL <- aln[order(aln$s_id),]


#################################
## Look at the data and format ##
#################################

## Let's figure out sites and stages data ...
alnagg <- aggregate(alnALL[("year")], alnALL[c("s_id", "flower","leaf","offset", "lat", "lon", "alt")],
                    FUN=length)




# PEP_ID seems unique

# Subset the data based on the above for now ...

alnuse20 <- aln


####################################
## Get mean at each site and plot ##
####################################

## summarizing data
if(mapsandsummaries){
  meanaln <-
    ddply(alnALL, c("s_id", "lon", "lat"), summarise,
          mean = mean(offset),
          mean.yr = mean(year),
          sd = sd(offset),
          sem = sd(offset)/sqrt(length(offset)))
}  
 
pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)

}


alnUSW <- pepnumber(aln, "s_id")

#################################################
## Fit hinge models for each species (in Stan) ##
#################################################


# now add hinge
aln<-alnUSW
aln$YEAR.hin <- aln$year
aln$YEAR.hin[which(aln$YEAR.hin<1980)] <- 1980
aln$YEAR.hin <- aln$YEAR.hin-1980
# Note that I tried centering and scaling doy and it did not speed up much. The 20 yrs model still said: 1000 transitions using 10 leapfrog steps per transition would take 215.36 seconds.


# Note to self: for betula lmer will fit random intercepts but not random slopes

N <- nrow(aln)
y <- aln$offset
J <- length(unique(aln$peporder))
sites <- aln$peporder
year <- aln$YEAR.hin
  # nVars <-1
  # Imat <- diag(1, nVars)
  
  fit.hinge.aln <- stan("stan/hinge_randslopesint.stan",
                        data=c("N","J","y","sites","year"), iter=2500, warmup=1500,
                        chains=4, cores=ncores)
  # control = list(adapt_delta = 0.95, max_treedepth = 15))

}  
save(fit.hinge.aln, file="stan/fit.hinge.aln.Rda")
sumera <- summary(fit.hinge.aln)$summary
sumera[grep("mu_", rownames(sumera)),] ###this give you the overall
mean(aln$year)
sumera<-as.data.frame(sumera)
sumera<-rownames_to_column(sumera, var = "peporder")



####plot
a <- 29.7411191
b <- 0.2780304
x0 <- 1960
x1 <- 1980

# graph
plot(c(1960,2010), c(0,69), type = "n", xlab = "", ylab = "", bty='l')
segments(x0=1960,y0=a, x1=1980,y1=a)
segments(x0=1980,y0=a, x1=2015,y1=a+b*35)





#########################
getstanpred <- function(dat, sitecolname, stansummary, predyear){
  siteslist <- unlist(unique(dat[sitecolname]))
  sumer.ints <- stansummary[grep("a\\[", rownames(stansummary)),]
  sumer.slopes <- stansummary[grep("b\\[", rownames(stansummary)),]
  stanfit <- data.frame(m=as.numeric(rep(NA, length(siteslist))),
                        pred=as.numeric(rep(NA, length(siteslist))), site=siteslist)
  for (sitehere in c(1:length(siteslist))){
    stanfit$m[sitehere] <- sumer.slopes[sitehere]
    stanfit$pred[sitehere] <- sumer.ints[sitehere]+sumer.slopes[sitehere]*predyear
  }
  return(stanfit)
}
predaln <- getstanpred(aln, "s_id", sumera, 3)

getlinpred <- function(dat, sitecolname, predyear){
  siteslist <- unique(dat[sitecolname])
  linfit <- data.frame(m=as.numeric(rep(NA, nrow(siteslist))),
                       pred=as.numeric(rep(NA, nrow(siteslist))))
  for (sitehere in c(1:nrow(siteslist))){
    subby <- subset(dat, s_id==siteslist[sitehere,])
    mod <- lm(offset~YEAR.hin, data=subby)
    linfit$m[sitehere] <- coef(mod)[2]
    linfit$pred[sitehere] <- coef(mod)[1]+coef(mod)[2]*predyear
  }
  return(linfit)
}

alnpred.lin <- getlinpred(aln, "s_id", 3)

plot(predaln$pred~alnpred.lin$pred, asp=1)
abline(lm(predaln$pred~alnpred.lin$pred))
ggplot(aln,aes(year,offset))+geom_point()+geom_abline(intercept=1960,slope=1)

#######just do it in brms
library(brms)
fit.aln.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aln) 
dat<-as.data.frame(coef(fit.aln.brms))
library(tibble)
dat<-rownames_to_column(dat, var = "peporder")
newdat<-merge(dat,aln)
colnames(newdat)
names(newdat)[1]<-"peporder"
names(newdat)[2]<-"Intercept"
names(newdat)[6]<-"slope"
newdata<-dplyr::select(newdat,peporder,Intercept,slope,s_id,year,offset,YEAR.hin)
newdata$hinge<-ifelse(newdata$YEAR.hin>0,1,0)
Enewdata$slope<-ifelse(newdata$YEAR.hin==0,0,newdata$slope)

alpha<-mean(newdata$Intercept)
beta<-mean(newdata$slope)
ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha,yend=alpha,x=1960, xend=1980),color="red")
ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha,yend=alpha+35*beta),color="red")

ggplot(newdata,aes(year,offset,group=peporder))+geom_segment(aes(y=Intercept,yend=Intercept,x=1960, xend=1980))+geom_segment(aes(y=alpha,yend=alpha,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=Intercept,yend=Intercept+35*slope))+geom_segment(aes(x=1980,xend=2015,y=alpha,yend=alpha+35*beta),color="red")+theme_base()+theme(legend.position="none")


plotty+geom_segment(aes(y=mean(Intercept),yend=alpha,x=1960, xend=1980),color="red")+geom_segment(aes(x=1980,xend=2015,y=mean(Intercept),yend=alpha+35*beta),color="red")
alpha<-mean(newdata$Intercept)
beta<-mean(newdata$slope)
alpha+35*beta




