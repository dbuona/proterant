### Explore the prunus data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library(tibble)
library(lubridate)

setwd("~/Documents/git/proterant/investment/input")
load("PrunusFLSs.Rda")
d1<-read.csv("midwest_round1Dec10.csv")
d<-read.csv("midwest_round1Dec11.csv")
##subset to flowering data
d.flo<-dplyr::filter(d,flowers=="Y")
d.flo<-dplyr::filter(d.flo,bbch.f<80)
d.flo<-dplyr::filter(d.flo,bbch.f>59)
table(d.flo$specificEpithet)

ggplot(d.flo,aes(bbch.v,bbch.f))+geom_point()+
  facet_wrap(~specificEpithet)+geom_smooth(method="lm")

####fuctions
extract_coefsF<-function(x){rownames_to_column(as.data.frame(coef(x, summary=TRUE,probs=c(0.025,0.25,0.75,0.975))),"species")
}


## basic model
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

### now control for date of observation as a covariate
d.flo<-filter(d.flo,!is.na(eventDate))
d.flo$eventDate2<-as.Date(d.flo$eventDate,format =  "%Y-%m-%d")
d.flo$doy<-yday(d.flo$eventDate2)

d.flo$doy.cent<-d.flo$doy-mean(d.flo$doy,na.rm=TRUE)

### run model
mod2<-brm(bbch.v~doy.cent+(doy.cent|specificEpithet),data=d.flo)
mod2.out<-extract_coefsF(mod2)
mod2.out<-dplyr::select(mod2.out,1:7)

colnames(mod2.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2.out$species <- factor(mod2.out$species, levels = mod2.out$species[order(mod2.out$Estimate)])

mod2a<-brm(bbch.v~doy.cent+(1|specificEpithet),data=d.flo)
mod2a.out<-extract_coefsF(mod2a)
mod2a.out<-dplyr::select(mod2a.out,1:7)

colnames(mod2a.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2a.out$species <- factor(mod2a.out$species, levels = mod2a.out$species[order(mod2a.out$Estimate)])

mod2 <- add_criterion(mod2, "loo",reloo=TRUE)
mod2a <- add_criterion(mod2a, "loo",reloo=TRUE)

loo_compare(mod2,mod2a) ### pooling on intercept model is better
jpeg("..//Plots/seeminglybestplot.jpeg")
ggplot(mod2a.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 19), breaks =  c(9,11,15,17,19))+ggtitle("Control for date")+
theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

pp_check(mod2a,nsamples=100)
d.flo$bbch.v.scale[d.flo$bbch.v==0]<-0
d.flo$bbch.v.scale[d.flo$bbch.v==7]<-1
d.flo$bbch.v.scale[d.flo$bbch.v==9]<-2
d.flo$bbch.v.scale[d.flo$bbch.v==10]<-3
d.flo$bbch.v.scale[d.flo$bbch.v==11]<-4
d.flo$bbch.v.scale[d.flo$bbch.v==14]<-5
d.flo$bbch.v.scale[d.flo$bbch.v==15]<-5
d.flo$bbch.v.scale[d.flo$bbch.v==17]<-6
d.flo$bbch.v.scale[d.flo$bbch.v==19]<-7

mod2a.scale<-brm(bbch.v.scale~doy.cent+(1|specificEpithet),data=d.flo)
pp_check(mod2a.scale)
mod2a.scale.out<-extract_coefsF(mod2a.scale)
mod2a.scale.out<-dplyr::select(mod2a.scale.out,1:7)

colnames(mod2a.scale.out)<-c("species","Estimate","SE","Q2.5","Q25","Q75","Q97.5")
mod2a.scale.out$species <- factor(mod2a.scale.out$species, levels = mod2a.scale.out$species[order(mod2a.scale.out$Estimate)])

jpeg("..//Plots/seeminglybestplot_rescaled.jpeg")
ggplot(mod2a.scale.out,aes(species,Estimate))+geom_point()+
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0,linetype="dotted")+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  scale_y_continuous(limits = c(0, 9), breaks =  c(2,4,5,6,7),labels=c(9,11,15,17,19))+ggtitle("Control for date, rescaled")+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))
dev.off()

hist(d.flo$bbch.v.scale)
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


?geom_rect()
save.image("PrunusFLSs.Rda")





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
