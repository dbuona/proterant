####FINAL PRUNUS ANAYSIS: see plummywork for scratch and hydraulic demand
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

graphics.off()
library(dplyr)
library(ggplot2)
library("rstan")
library(brms)


library(phytools)
library(ape)
library(lubridate)
library(stringr)
library("tidybayes")
library(raster)
library(sp)

#require(mapdata); require(maptools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/Documents/git/proterant/investment/Input")
load("TP_analyses_Rda")




palmer.b <- brick("..//Data/lbda-v2_kddm_pmdi_2017.nc")
#tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre") ##tree
d.flo<-read.csv("input_clean/FLS_clean.csv") ##data
tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre") ##tree

d.flo$species<-d.flo$specificEpithet ## whoops over wrote the id column here but we dont need if

#d.flo$logFLS<-log(d.flo$bbch.v.scale) ## make FLS linear

#d.flo$rtFLS<-sqrt(d.flo$bbch.v.scale) 
#ggplot(d.flo,aes(rtFLS))+geom_density()






d.flo.check <- d.flo %>%
  group_by(species) %>%
  filter(!(abs(doy - median(doy)) > 3*sd(doy)))

ggplot(d.flo,aes(doy))+geom_histogram(bins=20)+facet_wrap(~species)

ggplot(d.flo.check,aes(doy))+geom_histogram(bins=20)+facet_wrap(~species)
d.flo.check <- d.flo.check  %>% group_by(species) %>%  mutate(range = min(doy))
latpoints<-d.flo.check$lat #
lonpoints<-d.flo.check$lon

### Mean Temperature and precip
lonpoints<-lonpoints+360 # make vector of prunus coordinates
extract.pts<-cbind(lonpoints,latpoints)

extract.pts <- extract.pts[complete.cases(extract.pts), ]
extract.pts<-as.data.frame(extract.pts)
coordinates(extract.pts)<-~lonpoints+latpoints

p<-brick("~/Downloads/precip.mon.total.v501.nc")
t<-brick("~/Downloads/air.mon.mean.v501.nc")


mean.p <-calc(p, fun = mean,na.rm=TRUE) ## calculate the mean of each
mean.t <-calc(t, fun = mean,na.rm=TRUE)

Pp<-raster::extract(mean.p,extract.pts,method="simple")
Tt<-raster::extract(mean.t,extract.pts,method="simple")



Pp<-as.data.frame(Pp)
Tt<-as.data.frame(Tt)
joiner<-dplyr::select(d.flo.check,lat,lon,species)
joiner <- joiner[complete.cases(joiner),] 



rownames(Pp)<-rownames(joiner)
rownames(Tt)<-rownames(joiner)


goo<-cbind(joiner,Pp)
goo<-cbind(goo,Tt)

goo<-as.data.frame(goo)

palmer.b
palmer.b1 <-palmer.b[[1:1900]]## subset to pnly last century
palmer.b2 <-palmer.b[[1900:2018]]## subset to pnly last century
palmer.b1<-brick(palmer.b1)
palmer.b2<-brick(palmer.b2)

long <-calc(palmer.b1, fun = mean,na.rm=TRUE)
recent <-calc(palmer.b2, fun = mean,na.rm=TRUE)
lonpoints2<-d.flo.check$lon # make vector of prunus coordinates
latpoints2<-d.flo.check$lat
extract.pts2<-cbind(lonpoints2,latpoints2)

ext.recent<-raster::extract(recent,extract.pts2,method="simple")
ext.long<-raster::extract(long,extract.pts2,method="simple")


ext.recent<-as.data.frame(ext.recent)
ext.long<-as.data.frame(ext.long)
#ext.long<-filter(ext.long,!is.na(ext.long))

joiner2<-dplyr::select(d.flo.check,lat,lon,species)
#joiner2 <- joiner2[complete.cases(joiner2),] 

joiner2<-cbind(joiner2,ext.recent)

joiner2<-cbind(joiner2,ext.long)

PDSI <-joiner2 %>% group_by(species) %>% summarise(meanPDSIrec=mean(ext.recent,na.rm=TRUE),meanPDSIlong=mean(ext.long,na.rm=TRUE))
TP<-goo %>% group_by(species) %>% summarise(meanP=mean(Pp,na.rm=TRUE),meanT=mean(Tt,na.rm=TRUE))

compy<-left_join(PDSI,TP)

cor(compy$meanP,compy$meanPDSIrec)
cor(compy$meanP,compy$meanPDSIlong)

goober<-select(goo,species,lat,lon,Pp)
goober2<-select(joiner2,species,lat,lon,ext.recent,ext.long)
goober<-left_join(goober,goober2)
cor(goober$Pp,goober$ext.recent,use="pairwise.complete.obs")
cor(goober$Pp,goober$ext.long,use="pairwise.complete.obs")


jpeg("..//Plots/P_vs_pdsicorrs.jpeg",width=6,height=5,unit='in',res=250)
ggplot()+
  geom_point(data=compy,aes(meanPDSIlong,meanP),color="salmon")+geom_smooth(data=compy,aes(meanPDSIlong,meanP),method="lm",color="salmon",fullrange=TRUE,se=FALSE)+
  geom_point(data=compy,aes(meanPDSIrec,meanP),color="navy")+geom_smooth(data=compy,aes(meanPDSIrec,meanP),method="lm",color="navy",fullrange=TRUE,se=FALSE)+
  geom_text(data=compy,aes(meanPDSIrec,meanP,label=species),size=2.5,nudge_y = -.2,color="navy")+
geom_text(data=compy,aes(meanPDSIlong,meanP,label=species),size=2.5,nudge_y = .2,color="salmon")+ylab("pdsi")+
  xlab('precip')+ggthemes::theme_base()
dev.off()

    write.csv(TP,"input_clean/TP.csv")


####look at plasticity
extract.pts2<-as.data.frame(extract.pts)
varT<-raster::extract(t,extract.pts2,method="simple")
varT<-as.data.frame(varT)

varP<-raster::extract(p,extract.pts2,method="simple")
varP<-as.data.frame(varP)

joiney<-dplyr::select(d.flo.check,id,species,year,lon,lat)
joiney <- joiney[complete.cases(joiney),] 

vart<-cbind(joiney,varT)
vart<-tidyr::gather(vart,timeperiod,temperature,6:1421)

varp<-cbind(joiney,varP)
varp<-tidyr::gather(varp,timeperiod,precip,6:1421)



vart.t<-data.frame(timeperiod=unique(vart$timeperiod))


vart.t$year2<-rep(1900:2017,each=12)
vart.t$month<-rep(1:12,118)
vart<-left_join(vart,vart.t)
vart<-filter(vart,year==year2)
varty<- vart %>% group_by(id,year) %>% summarise(annualT=mean(temperature,na.rm=TRUE))



varp.p<-data.frame(timeperiod=unique(varp$timeperiod))
varp.p$year2<-rep(1900:2017,each=12)
varp.p$month<-rep(1:12,118)
varp<-left_join(varp,varp.p)
varp.summer<-filter(varp, month %in% c(6,7,8))
varpy.summer<- varp.summer %>% group_by(species) %>% summarise(annualP=mean(precip,na.rm=TRUE))
write.csv(varpy.summer,"Psummeronly.csv")
varp<-filter(varp,year==year2)
varpy<- varp %>% group_by(id,year) %>% summarise(annualP=mean(precip,na.rm=TRUE))



tempting<-left_join(d.flo.check,varty)
tempting<-left_join(tempting,varpy)


tree<-compute.brlen(tree, method = "Grafen")## make ultrametric
#check names
names.intree<-tree$tip.label # names the names
namelist<-unique(tempting$specificEpithet)
to.prune<-which(!names.intree%in%namelist) #prune the tree
pruned.by<-drop.tip(tree,to.prune)
plotTree(pruned.by)# this is the tree

###what are the tip labels in pruned phylogeny?

mytree.names<-pruned.by$tip.label # did i get them all
intersect(namelist,mytree.names) #yes

A <- ape::vcv.phylo(pruned.by) ## make acovarience matrix for brms models

cor(tempting$annualT,tempting$annualP,use = "pairwise.complete.obs")

save.image("..//Input/TP_analyses_Rda")

plasticT<-brm(bbch.v.scale~annualT+(annualT|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
plasticP<-brm(bbch.v.scale~annualP+(annualP|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

plasticT.doy<-brm(bbch.v.scale~doy+annualT+(annualT|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
plasticP.doy<-brm(bbch.v.scale~doy+annualP+(annualP|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

cor(tempting$doy,tempting$annualT,use="pairwise.complete.obs")
cor(tempting$doy,tempting$annualP,use="pairwise.complete.obs")
plasticT.doy2<-brm(bbch.v.scale~doy*annualT+(doy*annualT|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
plasticP.doy2<-brm(bbch.v.scale~doy*annualP+(doy*annualP|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))


coef(plasticP,probs = c(.055,.945,.25,.75))
fixef(plasticT,probs = c(.055,.945,.25,.75))
fixef(plasticT.doy,probs = c(.055,.945,.25,.75))

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
tempting$T.z<-zscore(tempting$annualT)
tempting$P.z<-zscore(tempting$annualP)
quantile(tempting$annualT,na.rm=TRUE,probs=(c(.15,.85)))

plasticTP<-brm(bbch.v.scale~T.z*P.z+(T.z*P.z|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
coef(plasticTP,probs =  c(.055,.945,.25,.75))
fixef(plasticT,probs =  c(.055,.945,.25,.75))
fixef(plasticP,probs =  c(.055,.945,.25,.75))
mean(tempting$annualT,na.rm=TRUE)

coef(plasticTP,probs =  c(.055,.945,.25,.75))
save.image("..//Input/TP_analyses_Rda")
new.data<-data.frame(species=rep(unique(tempting$species),each=2),annualT=rep(c(14.2,18.7),13))
new.data$specificEpithet<-new.data$species
fity<-fitted(plasticT,newdata = new.data,probs = c(.055,.945))
plotT<-cbind(new.data,fity)
prezest<- plotT %>%
  dplyr::select(annualT,species,contains("Estimate"))
colnames(prezest)[3:8]<-c(1,2,3,4,5,6)
prezest<-tidyr::gather(prezest,"stage","probability",3:8)

prezlow<- plotT %>%
  dplyr::select(annualT,species,contains("Q5.5"))
colnames(prezlow)[3:8]<-c(1,2,3,4,5,6)
prezlow<-tidyr::gather(prezlow,"stage","low",3:8)

prezhigh<- plotT %>%
  dplyr::select(annualT,species,contains("Q94.5"))
colnames(prezhigh)[3:8]<-c(1,2,3,4,5,6)
prezhigh<-tidyr::gather(prezhigh,"stage","high",3:8)

#prezhigh<-dplyr::select(prezhigh,-stage)
#prezlow<-dplyr::select(prezlow,-stage)

prezest<-left_join(prezest,prezlow)
prezest<-left_join(prezest,prezhigh)
prezest<-distinct(prezest)
prezest$stage<-as.numeric(prezest$stage)


prezest$bbch<-NA
prezest$bbch[which(prezest$stage==1 )]<- "BBCH 0"
prezest$bbch[which(prezest$stage==2 )]<- "BBCH 09"
prezest$bbch[which(prezest$stage==3)]<- "BBCH 11"
prezest$bbch[which(prezest$stage==4)]<- "BBCH 15"
prezest$bbch[which(prezest$stage==5 )]<- "BBCH 17"
prezest$bbch[which(prezest$stage==6 )]<- "BBCH 19"

ggplot(prezest,aes(bbch,probability))+
  geom_point(size=2,aes(color=as.factor(annualT),shape=as.factor(annualT)))+geom_ribbon(aes(x=as.numeric(stage),ymin=0,ymax=probability,fill=as.factor(annualT)),alpha=0.6)+facet_wrap(~species)+
  
#facet_wrap(~factor(species,levels=c("mexicana","umbellata","angustifolia","maritima","gracilis","munsoniana","alleghaniensis","nigra","americana","texana","rivularis","hortulana","subcordata")))+
  geom_errorbar(aes(ymin=low,ymax=high,color=as.factor(annualT)),width=0)+
  theme(strip.text = element_text(face = "italic"))+
ggthemes::theme_few()+
  scale_fill_viridis_d(option="C")+
  scale_color_manual(values=c("black","black"))+
   theme(strip.text = element_text(face = "italic"))+xlab("")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 270,vjust =.4,size = 7))



plotT<-select(plotT,-specificEpithet)

plotT<-tidyr::gather(plotT,"bbch","likelihood",3:5)

ggplot(plotT,aes(annualT,likelihood))+geom_point(aes(color=bbch))+facet_wrap(~species,scales="free_y")

fixef(plasticT)
coef(plasticT)



######
palmer.b1 <-palmer.b[[1:1900]]## subset to pnly last century
palmer.b2 <-palmer.b[[1900:2018]]## subset to pnly last century
palmer.b1<-brick(palmer.b1)
palmer.b2<-brick(palmer.b2)

long <-calc(palmer.b1, fun = mean,na.rm=TRUE)
recent <-calc(palmer.b2, fun = mean,na.rm=TRUE)
lonpoints2<-d.flo.check$lon # make vector of prunus coordinates
latpoints2<-d.flo.check$lat
extract.pts2<-cbind(lonpoints2,latpoints2)

#ext.recent<-raster::extract(recent,extract.pts2,method="simple")
#ext.long<-raster::extract(long,extract.pts2,method="simple")


x1b <- resample(long, recent)
dif <- recent - x1b

library("RColorBrewer")

par(mar = c(5, 5, 5, 5)) 
par(mfrow = c(1, 3))

jpeg("..//Plots/PDSIovertimemaps.jpeg",width=10,height=5,unit='in',res=300)
par(mar = c(5, 5, 5, 5)) 
par(mfrow = c(1, 3))
plot(long,col=rev( rainbow( 99, start=0,end=1 )),zlim=c(-.8,1.1))
plot(recent,  col=rev( rainbow( 99, start=0,end=1 ) ),zlim=c(-.8,1.1))
plot(dif,col=brewer.pal(n=11,name="PuOr"),zlim=c(-.6,1))
dev.off()    


### PDSI plasticity
#reset coordinate to this raster

ext$lat<-latpoints
ext$lon<-lonpoints2
colnames(ext)
pdsi.dater<-tidyr::gather(ext,"year","pdsi",1:119)
class(pdsi.dater$year)
pdsi.dater$year<-as.integer(pdsi.dater$year)

gner<-d.flo.check

pdsi.dater2<-dplyr::left_join(gner,pdsi.dater)
pdsi.dater2<-pdsi.dater2 %>% distinct()

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

pdsi.dater2$pdsi.z<-zscore(pdsi.dater2$pdsi)


plastic.PDSI<-brm(bbch.v.scale~pdsi+(pdsi|species)+(1|gr(specificEpithet, cov = A)),data=pdsi.dater2,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

fixef(plastic.PDSI,probs =  c(.055,.945,.25,.75))
fixef(plasticT,probs =  c(.055,.945,.25,.75))
fixef(plasticP,probs =  c(.055,.945,.25,.75))

PDSI.plas<-coef(plastic.PDSI,probs =  c(.055,.945,.25,.75))[1]
PDSI.plas<-as.data.frame(PDSI.plas)
PDSI.plas<-select(PDSI.plas,1:6)
colnames(PDSI.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
PDSI.plas$species<-rownames(PDSI.plas)

T.plas<-coef(plasticT,probs =  c(.055,.945,.25,.75))[1]
T.plas<-as.data.frame(T.plas)
T.plas<-dplyr::select(T.plas,1:6)
colnames(T.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
T.plas$species<-rownames(T.plas)


Tk<-fixef(plasticT,probs =  c(.055,.945,.25,.75))
T.plas<-as.data.frame(T.plas)
T.plas<-dplyr::select(T.plas,1:6)
colnames(T.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
T.plas$species<-rownames(T.plas)

P.plas<-coef(plasticP,probs =  c(.055,.945,.25,.75))[1]
P.plas<-as.data.frame(P.plas)
P.plas<-selecy(P.plas,1:6)
colnames(P.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
P.plas$species<-rownames(P.plas)

TP.plas<-coef(plasticTP,probs =  c(.055,.945,.25,.75))[1]
TP.plas<-as.data.frame(TP.plas)
T.plas<-dplyr::select(T.plas,1:6)
colnames(T.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
T.plas$species<-rownames(T.plas)

T.plas$sp<-factor(x=T.plas$species,levels=rev(c("mexicana","umbellata","angustifolia","maritima","gracilis","americana","munsoniana","alleghaniensis","nigra","hortulana","texana","rivularis","subcordata")))
T.plas<-left_join(T.plas,TP)

pp1<-ggplot(T.plas,aes(estimate,reorder(sp,meanT)))+geom_point(size=2.5)+geom_errorbarh(aes(xmin=Q5.5,xmax=Q94.5),size=.25,height=0)+
  geom_errorbarh(aes(xmin=Q25,xmax=Q75),height=0,size=1)+
  geom_vline(xintercept=0,linetype="dotdash")+xlab("temperature")+
  theme_minimal()+xlim(-.35,.4)+ylab("")+theme(axis.text.y = element_text(face="italic"))

pp2<-ggplot(P.plas,aes(estimate,reorder(species,-estimate)))+geom_point(size=2.5)+geom_errorbarh(aes(xmin=Q5.5,xmax=Q94.5),size=.25,height=0)+
  geom_errorbarh(aes(xmin=Q25,xmax=Q75),height=0,size=1)+geom_vline(xintercept=0,linetype="dotdash")+xlab("precipitation")+
  theme_minimal()+xlim(-.35,.4)+ylab("")+theme(axis.text.y = element_text(face="italic"))#+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())


pp3<-ggplot(PDSI.plas,aes(estimate,reorder(species,-estimate)))+geom_point(size=2.5)+geom_errorbarh(aes(xmin=Q5.5,xmax=Q94.5),size=.25,height=0)+
  geom_errorbarh(aes(xmin=Q25,xmax=Q75),height=0,size=1)+geom_vline(xintercept=0,linetype="dotdash")+xlab("PDSI")+
  theme_minimal()+xlim(-.35,.4)+ylab("")+theme(axis.text.y = element_text(face="italic"))#+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

pdf("..//Plots/pasticity_mus.pdf")
ggpubr::ggarrange(pp1,pp2,pp3,ncol=1,labels=(c("a)","b)","c)")))
dev.off()
