####plasticity
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

t<-brick("~/Downloads/air.mon.mean.v501.nc")

extract.pts2<-as.data.frame(extract.pts)
varT<-raster::extract(t,extract.pts2,method="simple")
varT<-as.data.frame(varT)

joiney<-dplyr::select(d.flo.check,id,species,year,lon,lat)
joiney <- joiney[complete.cases(joiney),] 

vart<-cbind(joiney,varT)
vart<-tidyr::gather(vart,timeperiod,temperature,6:1421)

vart.t<-data.frame(timeperiod=unique(vart$timeperiod))


vart.t$year2<-rep(1900:2017,each=12)
vart.t$month<-rep(1:12,118)
vart<-left_join(vart,vart.t)

vart<-filter(vart,year==year2)

vart.spring<-filter(vart, month %in% c(2,3,4))
varty<-vart.spring %>% group_by(id,year) %>% summarise(springT=mean(temperature,na.rm=TRUE))

varty2<-vart.spring %>%group_by(species) %>% summarise(meanspringT=mean(temperature,na.rm=TRUE))
write.csv(varty2,"input_clean/Tspring.csv")

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

plasticT<-brm(bbch.v.scale~springT+(springT|species)+(1|gr(specificEpithet, cov = A)),data=tempting,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

fixef(plasticT,probs = c(.055,.945,.25,.75))
coef(plasticT,probs = c(.055,.945,.25,.75))

T.plas<-coef(plasticT,probs =  c(.055,.945,.25,.75))[1]
T.plas<-as.data.frame(T.plas)
T.plas<-dplyr::select(T.plas,1:6)
colnames(T.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
T.plas$species<-rownames(T.plas)

T.plas$sp<-factor(x=T.plas$species,levels=rev(c("mexicana","umbellata","angustifolia","maritima","gracilis","americana","munsoniana","alleghaniensis","nigra","hortulana","texana","rivularis","subcordata")))

pp1<-ggplot(T.plas,aes(estimate,sp))+geom_point(size=2.5)+geom_errorbarh(aes(xmin=Q5.5,xmax=Q94.5),size=.5,height=0)+
  geom_errorbarh(aes(xmin=Q25,xmax=Q75),height=0,size=1)+
  geom_vline(xintercept=0,linetype="dotdash")+xlab(expression('Temperature'*~degree*C*''))+
  theme_minimal()+xlim(-.3,.1)+ylab("")+theme(axis.text.y = element_text(face="italic"))

save.image(file = "plasticity.Rda")


palmer.b2 <-palmer.b[[1900:2018]]#
palmer.b2<-brick(palmer.b2)

lonpoints2<-d.flo.check$lon # make vector of prunus coordinates
latpoints2<-d.flo.check$lat
extract.pts2<-cbind(lonpoints2,latpoints2)


ext<-raster::extract(palmer.b2,extract.pts2,method="simple")
colnames(ext)<-(c(1899:2017))
ext<-as.data.frame(ext)
ext$lat<-latpoints2
ext$lon<-lonpoints2
colnames(ext)
pdsi.dater<-tidyr::gather(ext,"year","pdsi",1:119)
class(pdsi.dater$year)
pdsi.dater$year<-as.integer(pdsi.dater$year)

gner<-d.flo.check

pdsi.dater2<-dplyr::left_join(gner,pdsi.dater)
pdsi.dater2<-pdsi.dater2 %>% distinct()





plastic.PDSI<-brm(bbch.v.scale~pdsi+(pdsi|species)+(1|gr(specificEpithet, cov = A)),data=pdsi.dater2,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

fixef(plastic.PDSI,probs =  c(.055,.945,.25,.75))
PDSI.plas<-coef(plastic.PDSI,probs =  c(.055,.945,.25,.75))[1]
PDSI.plas<-as.data.frame(PDSI.plas)
PDSI.plas<-dplyr::select(PDSI.plas,1:6)
colnames(PDSI.plas)<-c("estimate", "error","Q5.5","Q94.5","Q25","Q75")
PDSI.plas$species<-rownames(PDSI.plas)
PDSI.plas$sp<-factor(x=PDSI.plas$species,levels=rev(c("mexicana","umbellata","angustifolia","maritima","gracilis","americana","munsoniana","alleghaniensis","nigra","hortulana","texana","rivularis","subcordata")))

pp3<-ggplot(PDSI.plas,aes(estimate,sp))+geom_point(size=2.5)+geom_errorbarh(aes(xmin=Q5.5,xmax=Q94.5),size=.5,height=0)+
  geom_errorbarh(aes(xmin=Q25,xmax=Q75),height=0,size=1)+geom_vline(xintercept=0,linetype="dotdash")+xlab("PDSI")+
  theme_minimal()+xlim(-.3,.3)+ylab("")+theme(axis.text.y = element_text(face="italic"))#+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
jpeg("..//Plots/plastics.jpeg",width=5,height=7,unit='in',res=250)
ggpubr::ggarrange(pp1,pp3,ncol=1,labels = c("a)","b)"))
dev.off()
