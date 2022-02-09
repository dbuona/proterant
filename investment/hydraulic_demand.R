###water limation modeling and plotting
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

graphics.off()
library(dplyr)
library(ggplot2)
library(brms)
library("rstan")

library(phytools)
library(ape)
#library(lubridate)
library(stringr)
library("tidybayes")
library(raster)
library(bayesplot)

setwd("~/Documents/git/proterant/investment/Input")
load("hydro_demand.Rda")
hystscore<-read.csv("input_clean/hystscore.csv")
d<-read.csv("input_clean/pdsi_spi.csv")
d<-dplyr::select(d,-X)
d<-left_join(d,hystscore)

meanpdsi<-d %>% group_by(specificEpithet) %>% summarise(meanPDSI=mean(pdsi,na.rm=TRUE),sdPDSI=sd(pdsi,na.rm=TRUE),nPDSI=n(),sePDSI=sdPDSI / sqrt(nPDSI))
meanFLS<-read.csv("input_clean/FLSdescriptive.csv")
predy.bern<-read.csv("input_clean/bernoulliFLS_0_9.csv")

predy.bern.50<-dplyr::filter(predy.bern, quant=="50%")


meanpdsi<-filter(meanpdsi,specificEpithet %in% c(predy.bern.50$specificEpithet))
meany<-left_join(predy.bern.50,meanpdsi)

simp1<-brm(Estimate~meanPDSI,data=meany,family=Beta(link = "logit", link_phi = "log"))
summary(simp1)
simp1<-brm(Estimate|mi(Est.Error)~me(meanPDSI,sePDSI),data=meany,family=Beta(link = "logit", link_phi = "log"))

simp1<-brm(meanFLS~meanPDSI,data=meany)
simp2<-brm(meanFLS|mi(seFLS)~meanPDSI,data=meany)
simp3<-brm(meanFLS|mi(seFLS)~me(meanPDSI,sePDSI),warmup=3000,iter=4000,data=meany)

fixef(simp1,probs = c(.25,.75))
fixef(simp2,probs = c(.25,.75))
fixef(simp3,probs = c(.25,.75))

##simpleist model



###add phylogeny###
tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre") ##tree

##### give tree branch lengths
is.ultrametric(tree)
tree<-compute.brlen(tree, method = "Grafen")## make ultrametric
#check names
names.intree<-tree$tip.label # names the names
namelist<-unique(d$specificEpithet)
to.prune<-which(!names.intree%in%namelist) #prune the tree
pruned.by<-drop.tip(tree,to.prune)
plotTree(pruned.by)# this is the tree

###what are the tip labels in pruned phylogeny?

mytree.names<-pruned.by$tip.label # did i get them all
intersect(namelist,mytree.names) #yes

A <- ape::vcv.phylo(pruned.by)




mod.pdsi<-brm(pdsi~score+(1|specificEpithet),data=d,family=gaussian(),warmup=3000,iter=4000)
fixef(mod.pdsi,probs = c(0.25,0.75))

d$species<-d$specificEpithet

mod.pdsi.phylo<-brm(pdsi~score+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs

mod.pdsi.phylo.me<-brm(pdsi|mi(pdsi.sd)~score+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs


#mod.pdsisd.phylo<-brm(pdsi.sd~score+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs


get_variables(mod.pdsi.phylo)
lines<-mod.pdsi.phylo%>%
  spread_draws(b_Intercept,  b_score )

lines.nophylo<-mod.pdsi%>%
  spread_draws(b_Intercept,  b_score )

a<-ggplot()+
  geom_jitter(data=d,aes(score,pdsi),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=lines,aes(intercept=b_Intercept,slope=b_score),alpha=0.01,color="skyblue3")+
  geom_abline(data=lines,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="navy",size=2)+
  #geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean \nPDSI at \n collection sites")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)


### now covariate
d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")

d.fruit<-filter(d.fruit,fruit_type=="fleshy")

d.petal<-dplyr::select(d.petal,-X)
d.fruit<-dplyr::select(d.fruit,-X)

d.petal<-left_join(d.petal,hystscore)
d.fruit<-left_join(d.fruit,hystscore)

d.petal$species<-d.petal$specificEpithet
d.fruit$species<-d.fruit$specificEpithet


mod.petal<-brm(pental_lengh_mm~score+(1|id)+(1|specificEpithet),data=d.petal,family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.95)) 


grouppetal<-d.petal %>% group_by(id,specificEpithet,species) %>% summarise(meanind=mean(pental_lengh_mm,na.rm=TRUE))
grouppetal<-left_join(grouppetal,hystscore)


mod.petal.phylo<-brm(pental_lengh_mm~score+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) # 59 divergent transitions
mod.petalgr.phylo<-brm(meanind~score+(1|specificEpithet)+(1|gr(species, cov = A)),data=grouppetal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) # 8 divergent transitions


linespetal<-mod.petal.phylo%>%
  spread_draws(b_Intercept,  b_score )
linespetal.nophylo<-mod.petal%>%
  spread_draws(b_Intercept,  b_score )


mod.fruit<-brm(fruit_diam_mm~score+(1|id)+(1|specificEpithet),data=d.fruit,family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.95))
mod.fruit.phylo<-brm(fruit_diam_mm~score+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) # 5 diver

linesfruit<-mod.fruit.phylo%>%
  spread_draws(b_Intercept,  b_score )
linesfruit.nophylo<-mod.fruit%>%
  spread_draws(b_Intercept,  b_score )

b<-ggplot()+
  geom_jitter(data=d.fruit,aes(score,fruit_diam_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linesfruit,aes(intercept=b_Intercept,slope=b_score),alpha=0.01,color="skyblue3")+
  geom_abline(data=linesfruit,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="navy",size=2)+
  #geom_abline(data=linesfruit.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=linesfruit.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean \nFruit diameter")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)


c<-ggplot()+
  geom_jitter(data=d.petal,aes(score,pental_lengh_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linespetal,aes(intercept=b_Intercept,slope=b_score),alpha=0.01,color="skyblue3")+
  geom_abline(data=linespetal,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="navy",size=2)+
 # geom_abline(data=linespetal.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
#  geom_abline(data=linespetal.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean \nPetal Length")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)

e<-ggpubr::ggarrange(b,c,labels=c("b)","c)"))

jpeg("..//Plots/dataplots.jpeg", width=7, height=5,unit="in",res=300)
ggpubr::ggarrange(a,e,nrow=2,ncol=1,labels =c("a)" ))
dev.off()


####plasticity
d.um<-read.csv("input_clean/FLS_clean.csv") ##data
palmer.b <- brick("..//Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d.um$lon # make vector of prunus coordinates
latpoints<-d.um$lat #
extract.pts <- cbind(lonpoints,latpoints)
palmer.b
palmer.b2 <-palmer.b[[1900:2018]]## subset to pnly last century
palmer.b2<-brick(palmer.b2)
ext<-raster::extract(palmer.b2,extract.pts,method="simple")
colnames(ext)<-(c(1899:2017))
ext<-as.data.frame(ext)
ext$lat<-latpoints
ext$lon<-lonpoints
colnames(ext)
pdsi.dater<-tidyr::gather(ext,"year","pdsi",1:119)
class(pdsi.dater$year)

pdsi.dater$year<-as.integer(pdsi.dater$year)

head(pdsi.dater)
head(pdsi.dater)
joiner<-dplyr::select(d.um,specificEpithet,lat,year,lon,bbch.v.scale,doy)
head(joiner)
head(pdsi.dater)
pdsi.dater2<-dplyr::left_join(joiner,pdsi.dater)
pdsi.dater2<-pdsi.dater2 %>% distinct()
#brm(dry~(1|specificEpithet),data=pdsi.counter,family="bernoulli")

pdsi.dater2$pdsi.z<-zscore(pdsi.dater2$pdsi)
pdsi.dater2$doy.z<-zscore(pdsi.dater2$doy)
pdsi.dater2$species<-pdsi.dater2$specificEpithet
table(pdsi.dater2$species)
pdsi.dater2 %>%group_by(species) %>%summarize(mean(pdsi,na.rm=TRUE),sd(pdsi,na.rm=TRUE))

marty<-dplyr::filter(pdsi.dater2,species %in% c("alleghaniensis","angustifolia","nigra","texana","umbellata"))
table(pdsi.dater2$specificEpitet)
table(marty$species)
pdsi.dater2$species<-pdsi.dater2$specificEpithet

nrow(marty)
library(lme4)
goo<-brm(bbch.v.scale~doy+pdsi+(pdsi|species),data=marty)
goo2<-brm(bbch.v.scale~pdsi+doy+(pdsi|specificEpithet),data=marty)
coef(goo2,probs = c(.25,.75))

summary(goo)


varsp<-filter(pdsi.dater2,species %in% c("alleghaniensis","americana","angustifolia","gracilis","hortulana","maritima","munsoniana","nigra"))

##reprun


moda<-brm(bbch.v.scale~doy+pdsi+(pdsi|specificEpithet)+(1|gr(species, cov = A)),data=pdsi.dater2,data2=list(A=A),family=cumulative("logit"),control = list(adapt_delta=.99),warmup=3000,iter=4000)

summary(moda)

save.image("hydro_demand.Rda")

