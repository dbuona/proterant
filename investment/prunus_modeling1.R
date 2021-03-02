### Explore the prunus data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library(tibble)
library(lubridate)
library(stringr)
library("ncdf4")
library(raster)


setwd("~/Documents/git/proterant/investment/input")
load("PrunusFLSs.Rda")
d<-read.csv("midwest_round1Dec11.csv") # active datasheet

##subset to useful additions
d.add<-read.csv("species_additions - species_additions.csv")
d.add<-dplyr::filter(d.add,flowers=="Y")
unique(d.add$bbch.f)
##subset to flowering data
d.flo<-dplyr::filter(d,flowers=="Y")

###cleanup #clean up everybody
d.flo<-dplyr::filter(d.flo,bbch.f<80) #clean
d.flo<-dplyr::filter(d.flo,bbch.f>59) #clean
table(d.flo$specificEpithet)
table(d.add$specificEpithet)
d.flo<-rbind(d.flo,d.add)

d.flo$bbch.v<-ifelse(is.na(d.flo$bbch.v),0,d.flo$bbch.v)


###from all the stuff below we know a rescaled model with doy centered as a predictor is best
####fuctions
extract_coefsF<-function(x){rownames_to_column(as.data.frame(coef(x, summary=TRUE,probs=c(0.025,0.25,0.75,0.975))),"species")
}


### now control for date of observation as a covariate
#if month was missing we impute april
d.flo$month<-ifelse(d.flo$month==0,4,d.flo$month)
d.flo$month<-ifelse(is.na(d.flo$month),4,d.flo$month)

#if day was missing we imput 15
d.flo$day<-ifelse(d.flo$day==0,15,d.flo$day)
d.flo$day<-ifelse(is.na(d.flo$day),15,d.flo$day)

d.flo$year<-ifelse(d.flo$day>32,d.flo$day,d.flo$year)
d.flo$day<-ifelse(d.flo$day>32,15,d.flo$day)

#clean year and impute 1980 if missing
d.flo$year<-ifelse(d.flo$year==2,2002,d.flo$year)
d.flo$year<-ifelse(d.flo$year==5,2005,d.flo$year)
d.flo$year<-ifelse(d.flo$year==6,2006,d.flo$year)
d.flo$year<-ifelse(d.flo$year==198,1980,d.flo$year)
d.flo$year<-ifelse(is.na(d.flo$year),1980,d.flo$year)
d.flo$year<-ifelse(d.flo$year==0,1980,d.flo$year)
unique(d.flo$year)

d.flo$eventDate2<-paste(d.flo$year,d.flo$month,d.flo$day,sep="-")
d.flo<-filter(d.flo,!is.na(eventDate2))
d.flo$eventDate3<-as.Date(d.flo$eventDate2,format =  "%Y-%m-%d")
d.flo$doy<-yday(d.flo$eventDate3)

d.flo$doy.cent<-d.flo$doy-mean(d.flo$doy,na.rm=TRUE)
d.flo<-filter(d.flo,!is.na(doy.cent))
table(d.flo$specificEpithet)
## now rescale BBCH so its even and postive
unique(d.flo$bbch.v)
d.flo$bbch.v.scale<-NA
d.flo$bbch.v.scale[d.flo$bbch.v==0]<-1
d.flo$bbch.v.scale[d.flo$bbch.v==7]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==9]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==10]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==11]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==14]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==15]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==17]<-5
d.flo$bbch.v.scale[d.flo$bbch.v==19]<-6

table(d.flo$bbch.v.scale)

unique(d.flo$specificEpithet)

meanfls<-d.flo %>% dplyr::group_by(specificEpithet)%>% dplyr::summarise(meanFLS=mean(bbch.v.scale),sdFLS=sd(bbch.v.scale))

#d.flo$specificEpithet[d.flo$specificEpithet=="alleghaniensis"]<-"umbellata"
#d.flo$specificEpithet[d.flo$specificEpithet=="munsoniana"]<-"rivularis"






mod2a.scale<-brm(bbch.v.scale~doy.cent+(1|specificEpithet),data=d.flo)
pp_check(mod2a.scale,nsamples = 100)
 
mod2a.scale.out<-extract_coefsF(mod2a.scale)
mod2a.scale.out<-dplyr::select(mod2a.scale.out,1:7)

colnames(mod2a.scale.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2a.scale.out$species <- factor(mod2a.scale.out$species, levels = mod2a.scale.out$species[order(mod2a.scale.out$Estimate)])
write.csv(mod2a.scale.out,"fls_gaus_estimate.csv",row.names = FALSE)

jpeg("..//Plots/gaussplot_DAC.jpeg")
ggplot(mod2a.scale.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(1, 6),breaks=(1:6),labels = c(0,9,11,15,17,19))+
  ggtitle("Predicted vegetative stage during flowering")+ylab("BBCH")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()
d.flo$bbch.v.scale<-as.factor(d.flo$bbch.v.scale)

#d.flo<-dplyr::filter(d.flo,!is.na(bbch.v.scale))

d.flo$bbch.v.scale<-as.integer(d.flo$bbch.v.scale)

pdsi<-read.csv("raw_means_pdsi.csv")

colnames(pdsi)[1]<-"specificEpithet"
pdsi<-filter(pdsi, specificEpithet %in% unique(d.flo$specificEpithet))
d.twist<-left_join(d.flo,pdsi)


mod2a.ord<-brm(bbch.v.scale~doy.cent+me(meanpdsi,sdminpdsi)+(1|specificEpithet),data=d.twist,family=cumulative("logit"))
## predict it for average doy
pp_check(mod2a.ord,nsamples = 100)

new.data<-data.frame(specificEpithet=unique(d.flo$specificEpithet),doy.cent=rep(mean(d.flo$doy.cent),13))
predy<-fitted(mod2a.ord,newdata = new.data,probs = c(.025,.975))

predy<-as.data.frame(predy)
predy$species<-unique(d.flo$specificEpithet)


predy2<-predy %>%tidyr::gather("phase","likelihood",1:24)
predy.est<-filter(predy2,str_detect(phase, "^Estimate"))
predy.error<-filter(predy2,str_detect(phase, "^Q2.5"))
predy.error2<-filter(predy2,str_detect(phase, "^Q97.5"))

errorlow <- predy.error %>% 
  group_by(species) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species)

errorhigh <- predy.error2 %>% 
  group_by(species) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species)

result <- predy.est %>% 
  group_by(species) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species)

colnames(errorlow)[3]<-"Q25"
colnames(errorhigh)[3]<-"Q75"


result$Q25<-errorlow$Q25
result$Q75<-errorhigh$Q75
result1<-result

ggplot(data=result1,aes(species,likelihood))+geom_point()+facet_wrap(~phase)+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))




errorlow <- predy.error %>% 
  group_by(species,phase) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species)

errorhigh <- predy.error2 %>% 
  group_by(species,phase) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species)

result <- predy.est %>% 
  group_by(species,phase) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species)

colnames(errorlow)[3]<-"Q25"
colnames(errorhigh)[3]<-"Q75"


result$Q25<-errorlow$Q25
result$Q75<-errorhigh$Q75
result
result1$goo<-paste(result1$species,result1$phase)
result$goo<-paste(result$species,result$phase)

result$most<-NA
result$most<-ifelse(result$goo %in% c(result1$goo),"Y","N")

result$bbch<-NA
result$bbch[which(result$phase=="Estimate.1" )]<- "BBCH 0"
result$bbch[which(result$phase=="Estimate.2" )]<- "BBCH 09"
result$bbch[which(result$phase=="Estimate.3" )]<- "BBCH 11"
result$bbch[which(result$phase=="Estimate.4" )]<- "BBCH 15"
result$bbch[which(result$phase=="Estimate.5" )]<- "BBCH 17"
result$bbch[which(result$phase=="Estimate.6" )]<- "BBCH 19"

jpeg("..//Plots/orid_DAC.jpeg", width=11, height=4,unit="in",res=300)
ggplot(data=result,aes(species,likelihood))+geom_point(aes(color=most))+facet_wrap(~bbch,nrow=1)+
  geom_errorbar(aes(ymin=Q25,ymax=Q75,color=most),width=0)+ggthemes::theme_clean()+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
  
dev.off()


result$int<-NA
result$int[which(result$phase=="Estimate.1" )]<- 1
result$int[which(result$phase=="Estimate.2" )]<- 2
result$int[which(result$phase=="Estimate.3" )]<- 3
result$int[which(result$phase=="Estimate.4" )]<- 4
result$int[which(result$phase=="Estimate.5" )]<- 5
result$int[which(result$phase=="Estimate.6" )]<- 6


jpeg("..//Plots/orid_DAC2.jpeg", width=11, height=4,unit="in",res=300)
ggplot(data=result,aes(bbch,likelihood))+geom_point(aes(color=most))+
  geom_ribbon(aes(x=int,ymin=0,ymax=likelihood),fill="lightgray",alpha=0.6)+facet_wrap(~species,nrow=2)+
  geom_errorbar(aes(ymin=Q25,ymax=Q75,color=most),width=0)+ggthemes::theme_clean()+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

###################################
### Part II: cliamte change#################
########################################

d.peak<-filter(d.flo,bbch.f==65)


### hinge
d.peak<-transform(d.peak,YEAR.hin=ifelse(year<=1980,1980,year))

d.peak<-transform(d.peak,YEAR.hin=YEAR.hin-1980)

mod.flo.hing<-brm(doy~YEAR.hin+(YEAR.hin|specificEpithet),data=d.peak)

table(d$bbch.v)
d.bb<-filter(d,bbch.v %in% c(11,15))
d.bb<-transform(d.bb,YEAR.hin=ifelse(year<=1980,1980,year))
d.bb<-transform(d.bb,YEAR.hin=YEAR.hin-1980)
d.bb$eventDate2<-as.Date(d.bb$eventDate,format =  "%Y-%m-%d")
d.bb$doy<-yday(d.bb$eventDate2)


mod.clim.hing.bb<-brm(doy~YEAR.hin+(YEAR.hin|specificEpithet),data=d.bb)

fixef(mod.flo.hing)
fixef(mod.clim.hing.bb)
coef(mod.flo.hing)

new.data<-data.frame(specificEpithet=rep(unique(d.flo$specificEpithet),5),YEAR.hin=rep(c(0,10,20,30,40),each=13))
flopred<-predict(mod.flo.hing,newdata = new.data,probs = c(.025,.25,.75,.975))
flopred<-cbind(flopred,new.data)

leafpred<-predict(mod.clim.hing.bb,newdata = new.data,probs = c(.025,.25,.75,.975))
leafpred<-cbind(leafpred,new.data)

leafpred$phase<-"leafout"
flopred$phase<-"flowering"

predy<-rbind(leafpred,flopred)
pre<-filter(predy,YEAR.hin==0)
pre2<-pre
pre$year<-1900
pre2$year<-1980
pre<-rbind(pre,pre2)

predy$year<-NA
predy$year[which(predy$YEAR.hin==0)]<-1980
predy$year[which(predy$YEAR.hin==10)]<-1990
predy$year[which(predy$YEAR.hin==20)]<-2000
predy$year[which(predy$YEAR.hin==30)]<-2010
predy$year[which(predy$YEAR.hin==40)]<-2020
#predy<-rbind(predy,pre)
jpeg("..//Plots/climchange.jpeg",width = 16,height=8,units = 'in',res=200)
ggplot()+
  geom_linerange(data=pre,aes(y=Estimate,xmin=1900,xmax=1980,color=phase),size=1)+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Estimate,color=phase),size=1)+
  facet_wrap(~specificEpithet,nrow=2)+
  geom_smooth(method="lm",data=pre,aes(x=year,y=Q75,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+
  geom_smooth(method="lm",data=pre,aes(x=year,y=Q25,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+
  # geom_linerange(data=pre,aes(y=Q25,xmin=1900,xmax=1980,color=phase),size=.2,alpha=.2,linetype="dotted")+
  #geom_linerange(data=pre,aes(y=Q75,xmin=1900,xmax=1980,color=phase),size=.2,alpha=.2,linetype="dotted")+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Q25,color=phase),size=.4,alpha=.8,linetype="dashed",se = FALSE)+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Q75,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+  
  ggthemes::theme_base(base_size = 9)+scale_x_continuous(limits = c(1900, 2020), breaks =  c(1980,1980,2020))+
  scale_color_manual(values=c("darkorchid1","darkgreen"))
dev.off()

ggplot()+

  geom_smooth(method="lm",data=pre,aes(x=year,y=Estimate,color=phase,fill=phase),size=.4,alpha=.6,linetype="dashed",se=TRUE)+
  #geom_smooth(method="lm",data=pre,aes(x=year,y=Q75,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+
  #geom_smooth(method="lm",data=pre,aes(x=year,y=Q25,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+
  # geom_linerange(data=pre,aes(y=Q25,xmin=1900,xmax=1980,color=phase),size=.2,alpha=.2,linetype="dotted")+
  #geom_linerange(data=pre,aes(y=Q75,xmin=1900,xmax=1980,color=phase),size=.2,alpha=.2,linetype="dotted")+
  #geom_smooth(method="lm",data=predy,aes(x=year,y=Q25,color=phase),size=.4,alpha=.8,linetype="dashed",se = FALSE)+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Estimate,color=phase,fill=phase),size=.4,alpha=.6,linetype="dashed",se=TRUE)+  
  ggthemes::theme_base(base_size = 9)+scale_x_continuous(limits = c(1900, 2020), breaks =  c(1980,1980,2020))+
  scale_color_manual(values=c("darkorchid1","darkgreen"))+scale_fill_manual(values=c("darkorchid1","darkgreen"))


ggplot()+
  geom_smooth(method="lm",data=pre,aes(x=year,y=Estimate,color=phase),size=1,se=FALSE)+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Estimate,color=phase),size=1,se=FALSE)+
  geom_smooth(method="lm",data=pre,aes(x=year,y=Q75,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+
  geom_smooth(method="lm",data=pre,aes(x=year,y=Q25,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+
  # geom_linerange(data=pre,aes(y=Q25,xmin=1900,xmax=1980,color=phase),size=.2,alpha=.2,linetype="dotted")+
  #geom_linerange(data=pre,aes(y=Q75,xmin=1900,xmax=1980,color=phase),size=.2,alpha=.2,linetype="dotted")+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Q25,color=phase),size=.4,alpha=.8,linetype="dashed",se = FALSE)+
  geom_smooth(method="lm",data=predy,aes(x=year,y=Q75,color=phase),size=.4,alpha=.8,linetype="dashed",se=FALSE)+  
  ggthemes::theme_base(base_size = 9)+scale_x_continuous(limits = c(1900, 2020), breaks =  c(1980,1980,2020))+
  scale_color_manual(values=c("darkorchid1","darkgreen"))


##do we have enough data?
olddat<-read.csv("midwest_round1Dec10.csv")
d1<-dplyr::filter(olddat,flowers=="Y")
d1<-dplyr::filter(d1,bbch.f<80) #clean
d1<-dplyr::filter(d1,bbch.f>59) #clean
table(d1$specificEpithet)

### now control for date of observation as a covariate
d1<-filter(d1,!is.na(eventDate))
d1$eventDate2<-as.Date(d1$eventDate,format =  "%Y-%m-%d")
d1$doy<-yday(d1$eventDate2)

d1$doy.cent<-d1$doy-mean(d1$doy,na.rm=TRUE)
d1<-filter(d1,!is.na(doy.cent))

## now rescale BBCH so its even and postive
unique(d1$bbch.v)
d1$bbch.v.scale<-NA
d1$bbch.v.scale[d1$bbch.v==0]<-1
d1$bbch.v.scale[d1$bbch.v==7]<-2
d1$bbch.v.scale[d1$bbch.v==9]<-3
d1$bbch.v.scale[d1$bbch.v==10]<-4
d1$bbch.v.scale[d1$bbch.v==11]<-4
d1$bbch.v.scale[d1$bbch.v==14]<-5
d1$bbch.v.scale[d1$bbch.v==15]<-5
d1$bbch.v.scale[d1$bbch.v==17]<-6
d1$bbch.v.scale[d1$bbch.v==19]<-7

table(d1$bbch.v.scale)

mod2a.scale.loww<-brm(bbch.v.scale~doy.cent+(1|specificEpithet),data=d1)

pp_check(mod2a.scale)

mod2a.scale.out.low<-extract_coefsF(mod2a.scale.loww)
mod2a.scale.out.low<-dplyr::select(mod2a.scale.out.low,1:7)

colnames(mod2a.scale.out.low)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2a.scale.out.low$species <- factor(mod2a.scale.out.low$species, levels = mod2a.scale.out.low$species[order(mod2a.scale.out.low$Estimate)])
jpeg("..//Plots/gaussplot_DAC.low.jpeg")
ggplot(mod2a.scale.out.low,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(1, 7),breaks=(1:7),labels = c(0,7,9,11,15,17,19))+
  ggtitle("760 rows vs 926 for current")+ylab("BBCH")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

palmer.b <- brick("..//lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d$lon # make vector of prunus coordinates
latpoints<-d$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
sd.prunus <-calc(palmer.b, fun = sd,na.rm=TRUE) #
min.prunus <-calc(palmer.b, fun = min,na.rm=TRUE) 
ext<-raster::extract(mean.prunus,extract.pts,method="simple")
ext2<-raster::extract(sd.prunus,extract.pts,method="simple")
ext3<-raster::extract(min.prunus,extract.pts,method="simple")
d$pdsi<-ext
d$pdsi.sd<-ext2
d$pdsi.min<-ext3


dlist <- list(
  pdsi = d.flo$pdsi,
  pdsi.sd  = d.flo$pdsi.sdsd,
  species     = d.flo$specificEpithet
  ) 
 
mod.pdsi<-brm(data=d,pdsi~(1|specificEpithet))
mod.pdsi.min<-brm(data=d,pdsi.min~(1|specificEpithet))
#mod.pdsi.error<-brm(data=d,pdsi.min | mi(pdsi.sd)~(1|specificEpithet))


mod.min.out<-extract_coefsF(mod.pdsi.min)

colnames(mod.min.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod.min.out1<-filter(mod.min.out, species %in% unique(mod2a.scale.out$species))


mod.pdsi.out<-extract_coefsF(mod.pdsi)
summary(mod.pdsi)
colnames(mod.pdsi.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod.pdsi.out1<-filter(mod.pdsi.out, species %in% unique(mod2a.scale.out$species))
plot(mod.pdsi.out1$Estimate,mod2a.scale.out$Estimate)

jointgoo<-data.frame(FLS=mod2a.scale.out$Estimate, FLSlow=mod2a.scale.out$Q25,FLShigh=mod2a.scale.out$Q75,
           PDSI=mod.pdsi.out1$Estimate,PDSIlow=mod.pdsi.out1$Q25,PDSIhigh=mod.pdsi.out1$Q75,
           PDSImin=mod.min.out1$Estimate,PDSIminlow=mod.min.out1$Q25,PDSIminhigh=mod.min.out1$Q75)
jointgoo$species<-mod.pdsi.out1$species


## or should I just use the mins for each species
pdsi.means<-d %>%group_by((specificEpithet)) %>% summarise(meanpdsi=mean(pdsi,na.rm=TRUE),sdpdsi=sd(pdsi,na.rm=TRUE),minpdsi=mean(pdsi.min,na.rm=TRUE),sdminpdsi=sd(pdsi.min,na.rm=TRUE))


write.csv(jointgoo,"means_data.csv",row.names = FALSE)
write.csv(pdsi.means,"raw_means_pdsi.csv",row.names = FALSE)

jpeg("..//Plots/pdsireg.jpeg")
ggplot(jointgoo,aes(FLS,PDSI))+geom_point()+ylab("PDSI")+
  geom_errorbar(aes(x=FLS,ymin=PDSIlow,ymax=PDSIhigh))+
    scale_x_continuous(name = "FLS",1:5,breaks=c(1,2,3,4,5),labels = c("bbch 0","bbch 09","bbch 11","bbch 15","bbch 17"))+
  geom_errorbarh(aes(y=PDSI,xmin=FLSlow,xmax=FLShigh))+geom_smooth(method="lm",se = T,fullrange=T)+
 
  ggthemes::theme_few()
dev.off()

filter(jointgoo,species!="maritima")
jointgoo$species<-mod2a.scale.out$species
summary(lm(FLS~PDSI,data=jointgoo))

nooutlyer<-filter(jointgoo,species!="maritima")
ggplot(nooutlyer,aes(FLS,PDSI))+geom_point()+
  geom_errorbar(aes(x=FLS,ymin=PDSIlow,ymax=PDSIhigh))+
  geom_errorbarh(aes(y=PDSI,xmin=FLSlow,xmax=FLShigh))+geom_smooth(method="lm",se = F,fullrange=T)+
  geom_smooth(method="lm",aes(FLSlow,PDSIhigh),color="red",se = F,fullrange=T)+
  geom_smooth(method="lm",aes(FLShigh,PDSIlow),color="green",se = F,fullrange=T)+
  geom_smooth(method="lm",aes(FLShigh,PDSIhigh),color="pink",se = F,fullrange=T)+
  geom_smooth(method="lm",aes(FLSlow,PDSIlow),color="darkgreen",se = F,fullrange=T)+
  ggthemes::theme_few()
 
## or does justmin pdsi matter
minpdsi<- d.flo %>% group_by(specificEpithet) %>% summarise(minpdsi=min(pdsi,na.rm=TRUE))
colnames(minpdsi)[1]<-"species"
minpdsi<-left_join(minpdsi,jointgoo)

ggplot(minpdsi,aes(FLS,minpdsi))+geom_point()+geom_smooth(method="lm")
summary(lm(FLS~minpdsi,data=minpdsi))

ggplot(d.flo,aes(bbch.v.scale,pdsi))+geom_point()+geom_smooth(method="lm")
library(lme4)
summary(lmer(bbch.v.scale~doy.cent+pdsi+(1|specificEpithet),data=d.flo))

summary(lmer(bbch.v.scale~doy.cent+pdsi.extreme+(1|specificEpithet),data=d.flo))
save.image("PrunusFLSs.Rda")





##comapred partial and no pooling models
mod1<-brm(bbch.v~(1|specificEpithet),data=d.flo)
mod1.nopool<-brm(bbch.v~specificEpithet-1,data=d.flo)
summary(mod1.nopool)

newdata<-data.frame(specificEpithet=unique(d.flo$specificEpithet))
preda<-predict(mod1.nopool,newdata = newdata,probs = c(.025,.25,.75,.975))
newdata<-cbind(newdata,preda)
colnames(newdata)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
newdata$modtype<-"no pool"

mod1.out<-extract_coefsF(mod1)

colnames(mod1.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod1.out$species <- factor(mod1.out$species, levels = mod1.out$species[order(mod1.out$Estimate)])
mod1.out$modtype<-"partial pool"
mod1er.out<-rbind(mod1.out,newdata)
jpeg("..//Plots/poolingcomps.jpeg")
ggplot(mod1er.out)+geom_point(aes(species,Estimate,color=modtype))+
geom_errorbar(aes(x=species,ymin=Q2.5,ymax=Q97.5,color=modtype),width=0,linetype="dotted")+
  geom_errorbar(aes(x=species,ymin=Q25,ymax=Q75,color=modtype),width=0)+ggthemes::theme_clean()+
   scale_y_continuous(limits = c(-10, 30), breaks =  c(9,11,15,17,19))+ggtitle("Intercept only")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+facet_wrap(~modtype)
dev.off()
jpeg("../Plots/nocontrol.jpeg")
ggplot(mod1.out)+geom_point(aes(species,Estimate))+
  geom_errorbar(aes(x=species,ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(x=species,ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Intercept only")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()


### run model
##rerun for loo
mod1<-brm(bbch.v~(1|specificEpithet),data=d.flo)

mod2<-brm(bbch.v~doy.cent+(doy.cent|specificEpithet),data=d.flo)
mod2.out<-extract_coefsF(mod2)
mod2.out<-dplyr::select(mod2.out,1:7)

colnames(mod2.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2.out$species <- factor(mod2.out$species, levels = mod2.out$species[order(mod2.out$Estimate)])

mod2a<-brm(bbch.v~doy.cent+(1|specificEpithet),data=d.flo) ## this is ulitamte a good model
mod2a.out<-extract_coefsF(mod2a)
mod2a.out<-dplyr::select(mod2a.out,1:7)

colnames(mod2a.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2a.out$species <- factor(mod2a.out$species, levels = mod2a.out$species[order(mod2a.out$Estimate)])

mod1 <- add_criterion(mod1, "loo",reloo=TRUE)
mod2 <- add_criterion(mod2, "loo",reloo=TRUE)
mod2a <- add_criterion(mod2a, "loo",reloo=TRUE)

loo_compare(mod1,mod2,mod2a) ### pooling on intercept model is better
jpeg("..//Plots/seeminglybestplot.jpeg")
ggplot(mod2a.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Control for date")+
theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

pp_check(mod2a,nsamples=100)



###try a poison distrubution, not so much better
mod2a.pois<-brm(bbch.v.scale~doy.cent+(1|specificEpithet),data=d.flo,family=poisson())

pp_check(mod2a.pois)
mod2a.pois.out<-extract_coefsF(mod2a.pois)
colnames(mod2a.pois.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2a.pois.out2<-dplyr::select(mod2a.pois.out,2:7)
mod2a.pois.out2<-exp(mod2a.pois.out2)
mod2a.pois.out2$species<-mod2a.pois.out$species

mod2a.pois.out2$species <- factor(mod2a.pois.out2$species, levels = mod2a.pois.out2$species[order(mod2a.pois.out2$Estimate)])

jpeg("..//Plots/seeminglybestplot_poisson.jpeg")
ggplot(mod2a.pois.out2,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 7), breaks =  c(2,3,4,5,6,7),labels=c(9,10,11,15,17,19))+ggtitle("Control for date, poisson")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()







stop("not error, below is scrap")
d.flo$lat.cent<-d.flo$lat-mean(d.flo$lat,na.rm=TRUE)

### run model
mod3<-brm(bbch.v~doy.cent+lat.cent+(1|specificEpithet),data=d.flo)
mod3.out<-extract_coefsF(mod3)
mod3.out<-dplyr::select(mod3.out,1:7)

colnames(mod3.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod3.out$species <- factor(mod3.out$species, levels = mod3.out$species[order(mod3.out$Estimate)])


c<-ggplot(mod3.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Control for date and lat")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))

ggpubr::ggarrange(a,b,c,nrow=2,ncol=2)

#### let day and lat vary by sp
mod4<-brm(bbch.v~doy.cent+lat.cent+(doy.cent+lat.cent|specificEpithet),data=d.flo)
mod4.out<-extract_coefsF(mod4)
mod4.out<-dplyr::select(mod4.out,1:7)

colnames(mod4.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod4.out$species <- factor(mod4.out$species, levels = mod4.out$species[order(mod4.out$Estimate)])

dd<-ggplot(mod4.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Control for date,lat, random slope")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))

ggpubr::ggarrange(a,b,c,d,nrow=2,ncol=2)

mod5<-brm(bbch.v~doy.cent*lat.cent+(doy.cent*lat.cent|specificEpithet),data=d.flo)
mod5.out<-extract_coefsF(mod5)
mod5.out<-dplyr::select(mod5.out,1:7)

colnames(mod5.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod5.out$species <- factor(mod5.out$species, levels = mod5.out$species[order(mod5.out$Estimate)])

e<-ggplot(mod5.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Control for date x lat, random slope")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))



mod1 <- add_criterion(mod1, "loo")
mod2 <- add_criterion(mod2, "loo")
mod3 <- add_criterion(mod3, "loo")
mod4 <- add_criterion(mod4, "loo")
mod5 <- add_criterion(mod5, "loo")
library(loo)
loo1<-loo(mod1) #-2479.4
loo2<-loo(mod2) # -2361.7 
loo3<-loo(mod3) #-2362.8 
loo4<-loo(mod4) #-2364.1
loo5<-loo(mod5) #-2364.6

loo_compare(loo3, loo2,loo4,loo5)
##model 2 is best for loo
?loo_compare()
pp_check(mod1) ### this might be why we need to rescale values so they are even

###try recoded model
unique(d.flo$bbch.v)
d$.flo$bbch.v
d.flo$bbch.v.scale[d.flo$bbch.v==0]<-0
d.flo$bbch.v.scale[d.flo$bbch.v==7]<-1
d.flo$bbch.v.scale[d.flo$bbch.v==9]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==10]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==11]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==14]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==15]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==17]<-5
d.flo$bbch.v.scale[d.flo$bbch.v==19]<-6

### try mod 1 again
mod1a<-brm(bbch.v.scale~(1|specificEpithet),data=d.flo)
pp_check(mod1a)

mod1a.out<-extract_coefsF(mod1a)
colnames(mod1a.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod1a.out$species <- factor(mod1a.out$species, levels = mod1a.out$species[order(mod1a.out$Estimate)])



aa<-ggplot(mod1a.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 6))+ggtitle("Intercept only")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
ggpubr::ggarrange(a,aa)

mod1a <- add_criterion(mod1a, "loo")



loo_compare(mod1,mod1a)

##compare to pdsi estiamtes ###No longer scrap. Need to arrange
palmer.b <- brick("..//lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d.flo$lon # make vector of prunus coordinates
latpoints<-d.flo$lat #
extract.pts <- cbind(lonpoints,latpoints)

mean.prunus <-calc(palmer.b, fun = mean,na.rm=TRUE) #average palmer drought index acrosss time
ext<-raster::extract(mean.prunus,extract.pts,method="simple")

d.flo$pdsi<-ext

mod.pdsi<-brm(pdsi~(1|specificEpithet),data=d.flo)

mod.pdsi.out<-extract_coefsF(mod.pdsi)

colnames(mod.pdsi.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")

jpeg("..//Plots/pdsi.jpeg",width = 10,height=8,units = 'in',res=200)
ggplot()+
  geom_point(data=mod2a.out,aes(species,Estimate),color="black")+
  geom_errorbar(data=mod2a.out,aes(x=species,ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted",color="black")+
  geom_errorbar(data=mod2a.out,aes(x=species,ymin=Q25,ymax=Q75),width=0,color="black")+ggthemes::theme_clean()+
  scale_y_continuous( breaks =  c(9,11,15,17,19), sec.axis=sec_axis(~.-12, name="PDSI"))+
  geom_point(data=mod.pdsi.out,aes(species,Estimate*10+12),color="firebrick")+
  geom_errorbar(data=mod.pdsi.out,aes(x=species,ymin=Q2.5*10+12,ymax=Q97.5*10+12),width=0,linetype="dotted",color="firebrick")+
  geom_errorbar(data=mod.pdsi.out,aes(x=species,ymin=Q25*10+12,ymax=Q75*10+12),width=0,color="firebrick")+ggthemes::theme_clean()+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

d.flotax<-d.flo

d.flotax$specificEpithet<-ifelse(d.flotax$specificEpithet=="munsoniana","rivularis",d.flotax$specificEpithet)
unique(d.flotax$specificEpithet)
d.flotax$specificEpithet<-ifelse(d.flotax$specificEpithet=="alleghaniensis","umbellata",d.flotax$specificEpithet)

mod2atax<-brm(bbch.v~doy.cent+(1|specificEpithet),data=d.flotax)
mod2atax.out<-extract_coefsF(mod2atax)
mod2atax.out<-dplyr::select(mod2atax.out,1:7)

colnames(mod2atax.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2atax.out$species <- factor(mod2atax.out$species, levels = mod2atax.out$species[order(mod2atax.out$Estimate)])




jpeg("..//Plots/seeminglybestplottaxa.jpeg")
ggplot(mod2atax.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Control for date, limited taxa")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

### oridianl https://kevinstadler.github.io/blog/bayesian-ordinal-regression-with-random-effects-using-brms/
d.flo$bbch.v.scale<-as.numeric(d.flo$bbch.v.scale)
d.flo$bbch.v.scale<-d.flo$bbch.v.scale+1
unique(d.flo$bbch.v.scale)



#1=0
#2=7
#3=9
#4=11
#5=15
#6=17
#7=19

ppa<-pp_check(mod2a.scale,nsamples = 100)
ppb<-pp_check(mod2a.ord,nsamples=100)

jpeg("..//Plots/ordinal_reg.jpeg")
ggpubr::ggarrange(ee,cc,ppa,ppb,nrow=2,ncol=2,heights = c(0.7,0.3))
dev.off()




