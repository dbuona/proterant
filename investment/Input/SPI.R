## extrac spi Processing of the Subset Data from RDA dataset ds298.0 - 'Standardized Precipitation Index (SPI) for Global Land Surface (1949-2012)
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


setwd("~/Documents/git/proterant/investment/")
spi <- brick("Data/470067.spi12.spi3_6_12_1deg_cru_ts_3_21_1949_2012.nc")
spi




###for spi
indices <- format(as.Date(names(spi), format = "X%Y.%m.%d"), format = "%m")
indices <- as.numeric(indices)
n<-names(spi)
nn <- as.integer(substr(n,7,8))

subj <-raster::subset(spi, which(nn %in% c(1,2,3,4,5,6)))
subj2 <-raster::subset(spi, which(nn %in% c(7,8,9,10,11,12)))
names(subj)
Monthspi<- stackApply(spi, indices, fun = mean)

minyear <-calc(spi, fun = min,na.rm=TRUE) 
minspring <-calc(subj, fun = min,na.rm=TRUE) 
minsummer <-calc(subj2, fun = min,na.rm=TRUE) 

meanyear <-calc(spi, fun = mean,na.rm=TRUE) 
meanspring <-calc(subj, fun = mean,na.rm=TRUE)
meansummer <-calc(subj2, fun = mean,na.rm=TRUE)


d<-read.csv("Input/input_clean/pruno_clean_pdsi.csv")

lonpoints<-d$lon # make vector of prunus coordinates
latpoints<-d$lat #



extract.pts <- cbind(lonpoints,latpoints)

ext<-raster::extract(minyear,extract.pts,method="simple")
ext2<-raster::extract(minspring,extract.pts,method="simple")
ext2a<-raster::extract(minsummer,extract.pts,method="simple")

ext3<-raster::extract(meanyear,extract.pts,method="simple")
ext4<-raster::extract(meanspring,extract.pts,method="simple")
ext4a<-raster::extract(meanspring,extract.pts,method="simple")


d$spiyear.min<-ext
d$spispring.min<-ext2
d$spisummer.min<-ext2a

d$spiyear.mean<-ext3
d$spispring.mean<-ext4
d$spisummer.mean<-ext4a


write.csv(d,"Input/input_clean/pdsi_spi.csv")

spi.means<-d %>%dplyr::group_by(specificEpithet) %>% dplyr::summarise(minspi=mean(spiyear.min,na.rm=TRUE),sdspi=sd(spiyear.min,na.rm=TRUE))
spispring.means<-d %>%group_by(specificEpithet) %>% summarise(minspi=mean(spispring.min,na.rm=TRUE),sdspi=sd(spispring.min,na.rm=TRUE))

d<-dplyr::left_join(d,spi.means)
d$withinspi<-d$spiyear.min-d$minspi
d$logFLS<-log(d$bbch.v.scale)


###add phylogeny

tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre")
is.ultrametric(tree)
tree<-compute.brlen(tree, method = "Grafen")
## make ultrametric

names.intree<-tree$tip.label # names the names
namelist<-unique(d$specificEpithet)

to.prune<-which(!names.intree%in%namelist) #prun the tree
pruned.by<-drop.tip(tree,to.prune)
plotTree(pruned.by)# this is the tree


###what are the tip labels in pruned phylogeny?

mytree.names<-pruned.by$tip.label # did i get them all
intersect(namelist,mytree.names) #yes



A <- ape::vcv.phylo(pruned.by)
plot(pruned.by)

##zscore


d$spiyear.min.z<-zscore(d$spiyear.min)
d$spiyear.mean.z<-zscore(d$spiyear.mean)
d$spispring.min.z<-zscore(d$spispring.min)
d$spispring.mean.z<-zscore(d$spispring.mean)

spi.mod.z<-brm(spiyear.min.z~(1|specificEpithet),data=d,warmup=2500,iter=4000)
pp_check(spi.mod.z)
spi.mod.z.phylo<-brm(spiyear.min.z~(1|gr(specificEpithet, cov = A)),data=d,data2=list(A=A),warmup=3000,iter=4000)
coef(spi.mod.z)
coef(spi.mod.z.phylo)
#petalmod.z<- brm(petal.z~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
petalmod.z.phylo<- brm(petal.z~(1|id)+(1|gr(specificEpithet, cov = A)),data=d.petal,data2=list(A=A),warmup=2500,iter=4000)

coef(petalmod.z)
coef(petalmod.z.phylo)
fruitmod.z<- brm(fruit.z~(1|id)+(1|specificEpithet),data=d.fruit,warmup=2500,iter=4000)
fruitmod.z.phylo<- brm(fruit.z~(1|id)+(1|gr(specificEpithet, cov = A)),data=d.fruit,data2=list(A=A),warmup=2500,iter=4000)
coef(fruitmod.z)
coef(fruitmod.z.phylo)


#phenmod.z<- brm(phen.z~(1|specificEpithet),data=d.phen,warmup=2500,iter=4000)

pdsiout<-dplyr::select(as.data.frame(coef(pdsi.mod.z.phylo)),1:2)
petalout<-dplyr::select(as.data.frame(coef(petalmod.z.phylo)$specificEpithet),1:2)
fruitout<-dplyr::select(as.data.frame(coef(fruitmod.z.phylo)$specificEpithet),1:2)


#spi.mins<-d %>%group_by((specificEpithet)) %>% summarise(minspi=min(spiyear,na.rm=TRUE),sdspi=sd(spiyear,na.rm=TRUE))
#spispring.mins<-d %>%group_by((specificEpithet)) %>% summarise(minspi=min(spispring,na.rm=TRUE),sdspi=sd(spispring,na.rm=TRUE))


meanpetals<-d.petal %>% group_by(specificEpithet) %>% summarise(meanPETAL=mean(pental_lengh_mm,na.rm=TRUE),sdPETAl=sd(pental_lengh_mm,na.rm=TRUE))
meanfruit<-d.fruit %>% group_by(specificEpithet) %>% summarise(meanFRUIT=mean(fruit_diam_mm,na.rm=TRUE),sdFRUIT=sd(fruit_diam_mm,na.rm=TRUE))

d<-left_join(d,meanpetals)
d<-left_join(d,meanfruit)



d$petal.z<-zscore(d$meanPETAL)
d$fruit.z<-zscore(d$meanFRUIT)
d$doy.z<-zscore(d$doy)

##also do it to the respone for later option
d$FLS.z<-zscore(d$bbch.v.scale)
d$lat.z<-zscore(d$lat)



library(lme4)
mod1<-lmer(logFLS~petal.z+fruit.z+spiyear.min.z+doy.z+(doy.z|specificEpithet),data=d)
mod1a<-lmer(logFLS~petal.z+fruit.z+spispring.min.z+doy.z+(doy.z|specificEpithet),data=d)
mod1b<-lmer(logFLS~petal.z+fruit.z+spiyear.mean.z+doy.z+(doy.z|specificEpithet),data=d)
mod1c<-lmer(logFLS~petal.z+fruit.z+spispring.mean.z+doy.z+(doy.z|specificEpithet),data=d)

summary(mod1)
summary(mod1a)
summary(mod1b)
summary(mod1c)
mod1.ord<-brm(bbch.v.scale~doy.z+petal.z+fruit.z+spiyear.min.z+(doy.z|specificEpithet),data=d,family=cumulative("logit"))
mod1a.ord<-brm(bbch.v.scale~doy.z+petal.z+fruit.z+spispring.min.z+(doy.z|specificEpithet),data=d,family=cumulative("logit"))
mod1.ord.ranef<-brm(bbch.v.scale~doy.z+petal.z+fruit.z+spiyear.min.z+(doy.z+spiyear.min.z|specificEpithet),data=d,family=cumulative("logit"))

mod1.ord<-add_criterion(mod1.ord,"")
?add_criterion()
mod1a.ord
mod1.ord.ranef

loo_compare(mod1.ord,mod1a.ord,mod1.ord.ranef)

library(tidybayes)
output<-mod1.ord %>%
  spread_draws(b_petal.z,b_fruit.z,b_spiyear.min.z)
colnames(output)
output <-output %>% tidyr::gather("var","estimate",4:6)
library(bayesplot)

ggplot(output,aes(y = var, x = estimate)) +
  tidybayes::stat_halfeye(fill="darkorchid1")+ggthemes::theme_few()+
  geom_vline(xintercept=0,linetype="dashed")



fixef(mod1.ord,probs = c(0.25,.75))
fixef(mod1a.ord,probs = c(0.25,.75))
fixef(mod1.ord.ranef,probs = c(0.25,.75))
coef(mod1.ord.ranef,probs = c(0.25,.75))








mod.drought<-lmer(spiyear.min~FLS.z+doy.z+(doy.z+FLS.z|specificEpithet),data=d)
summary(mod.drought)

mod.drought2<-lmer(spispring.min~doy.z+FLS.z+petal.z+fruit.z+(1|specificEpithet),data=d)
mod.drought3<-lmer(spiyear.mean~doy.z+FLS.z+petal.z+fruit.z+(1|specificEpithet),data=d)
mod.drought4<-lmer(spispring.mean~doy.z+FLS.z+petal.z+fruit.z+(1|specificEpithet),data=d)

summary(mod.drought2)
summary(mod.drought3)
summary(mod.drought4)

mod.bayes.drought<-brm(spiyear.min~doy.z+FLS.z+me(petal.z+fruit.z+lat.z+(FLS.z+petal.z+fruit.z|specificEpithet),data=d)
fixef(mod.bayes.drought,probs = c(0.25,.75))
coef(mod.bayes.drought)
