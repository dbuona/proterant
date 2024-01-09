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
lonpoints<-lonpoints+360 # make vector of prunus coordinates



extract.pts<-cbind(lonpoints,latpoints)

extract.pts <- extract.pts[complete.cases(extract.pts), ]
extract.pts<-as.data.frame(extract.pts)
coordinates(extract.pts)<-~lonpoints+latpoints

p<-brick("~/Downloads/precip.mon.total.v501.nc")
t<-brick("~/Downloads/air.mon.mean.v501.nc")


mean.p <-calc(p, fun = mean,na.rm=TRUE)
mean.t <-calc(t, fun = mean,na.rm=TRUE)

Pp<-raster::extract(mean.p,extract.pts,method="simple")
Tt<-raster::extract(mean.t,extract.pts,method="simple")
#cellStats(climy,stat = "mean")
#names(climy) <- c("Temp","Prec")
#tp<-extract(climy,extract.pts,method="simple")
Pp<-as.data.frame(Pp)
Tt<-as.data.frame(Tt)
joiner<-dplyr::select(d.flo.check,lat,lon,species)
joiner <- joiner[complete.cases(joiner),] 



rownames(Pp)<-rownames(joiner)
rownames(Tt)<-rownames(joiner)


goo<-cbind(joiner,Pp)
goo<-cbind(goo,Tt)

goo<-as.data.frame(goo)




TP<-goo %>% group_by(species) %>% summarise(meanP=mean(Pp,na.rm=TRUE),meanT=mean(Tt,na.rm=TRUE))

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
#vart<-arrange(vart,timeperiod)
#vart.t<-data.frame(timeperiod=vart$timeperiod)


vart.t$year2<-rep(1900:2017,each=12)
vart.t$month<-rep(1:12,118)
vart<-left_join(vart,vart.t)
vart<-filter(vart,year==year2)
varty<- vart %>% group_by(id,year) %>% summarise(annualT=mean(temperature,na.rm=TRUE))

tempting<-left_join(d.flo.check,varty)

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
plasticT<-brm(bbch.v.scale~annualT+(annualT|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
quantile(tempting$annualT,na.rm=TRUE,probs=(c(.15,.85)))

new.data<-data.frame(species=rep(unique(tempting$species),each=2),annualT=rep(c(9,18),13))
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










palmer.b
palmer.b1 <-palmer.b[[1:1900]]## subset to pnly last century
palmer.b2 <-palmer.b[[1900:2018]]## subset to pnly last century
palmer.b1<-brick(palmer.b1)
palmer.b2<-brick(palmer.b2)

long <-calc(palmer.b1, fun = mean,na.rm=TRUE)
recent <-calc(palmer.b2, fun = mean,na.rm=TRUE)

x1b <- resample(long, recent)
dif <- recent - x1b

library("RColorBrewer")

par(mar = c(5, 5, 5, 5)) 
par(mfrow = c(1, 3))

plot(long,col=rev( rainbow( 99, start=0,end=1 )),zlim=c(-.8,1.1))
plot(recent,  col=rev( rainbow( 99, start=0,end=1 ) ),zlim=c(-.8,1.1))
plot(dif,col=brewer.pal(n=11,name="PuOr"),zlim=c(-.6,1))
    


ext<-raster::extract(palmer.b2,extract.pts,method="simple")
colnames(ext)<-(c(1899:2017))
ext<-as.data.frame(ext)
ext$lat<-latpoints
ext$lon<-lonpoints
colnames(ext)
pdsi.dater<-tidyr::gather(ext,"year","pdsi",1:119)
class(pdsi.dater$year)
pdsi.dater$year<-as.integer(pdsi.dater$year)

joiner<-d.flo.check

pdsi.dater2<-dplyr::left_join(joiner,pdsi.dater)
pdsi.dater2<-pdsi.dater2 %>% distinct()

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

pdsi.dater2$pdsi.z<-zscore(pdsi.dater2$pdsi)
library(geodata)
#climy<-worldclim_country("United States",lon=lonpoints, lat=latpoints,var="bio",res=0.5,path="~/Documents/git/proterant/investment/Input/")
#climy$wc2.1_30s_bio_17
climy <- getData("worldclim",lon=lonpoints, lat=latpoints ,var="bio",res=2.5)
climy <- climy[[c(12)]]
climy<-calc(climy, fun = mean,na.rm=TRUE)
geodata_path()






crs(mean.p)<-crs(palmer.b)




#wa#pdsi.dater2$doy.z<-zscore(pdsi.dater2$doy)
#pdsi.dater2$species<-pdsi.dater2$specificEpithet
#table(pdsi.dater2$species)
#pdsi.dater2 %>%group_by(species) %>%summarize(mean(pdsi,na.rm=TRUE),sd(pdsi,na.rm=TRUE))

plastic.mod.ord.4review.nooutlier2<-brm(bbch.v.scale~pdsi+(pdsi|species)+(1|gr(specificEpithet, cov = A)),data=pdsi.dater2,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
coef(plastic.mod.ord.4review.nooutlier2,probs =c(0.055,.945))
table(pdsi.dater2$specificEpithet)
nigra<-filter(pdsi.dater2,species=="nigra")

plastic.mod.ngra<-brm(bbch.v.scale~pdsi,data=nigra,family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

#plastic.mod<-brm(bbch.v.scale~doy.z+pdsi.z+(pdsi.z|specificEpithet)+(1|gr(species, cov = A)),data=pdsi.dater2,data2=list(A=A),family=cumulative("logit"),control = list(adapt_delta=.99),warmup=3000,iter=4000)

#plastic.mod.noslp<-brm(bbch.v.scale~doy.z+pdsi.z+(1|specificEpithet)+(1|gr(species, cov = A)),data=pdsi.dater2,data2=list(A=A),family=cumulative("logit"),control = list(adapt_delta=.99),warmup=3000,iter=4000)


fixef(plastic.mod,probs = c(.25,.75))
fixef(plastic.mod.noslp,probs = c(.25,.75))




#air<-brick("air.mon.ltm.v501.nc")
air2<-brick("air.mon.mean.v501.nc")
#crs(palmer.b2) <- crs(air2)
#crs(air2) <- crs(palmer.b2)
air2
plot(air2)
air.extract<-extract(air2,extract.pts)

