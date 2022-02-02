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
hystscore<-read.csv("input_clean/hystscore.csv")
d<-read.csv("input_clean/pdsi_spi.csv")
d<-dplyr::select(d,-X)
d<-left_join(d,hystscore)

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



get_variables(mod.pdsi.phylo)
lines<-mod.pdsi.phylo%>%
  spread_draws(b_Intercept,  b_score )

lines.nophylo<-mod.pdsi%>%
  spread_draws(b_Intercept,  b_score )

a<-ggplot()+
  geom_jitter(data=d,aes(score,pdsi),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=lines,aes(intercept=b_Intercept,slope=b_score),alpha=0.01,color="navy")+
  geom_abline(data=lines,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="navy",size=2)+
  geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
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
  geom_abline(data=linesfruit,aes(intercept=b_Intercept,slope=b_score),alpha=0.01,color="navy")+
  geom_abline(data=linesfruit,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="navy",size=2)+
  geom_abline(data=linesfruit.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  geom_abline(data=linesfruit.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean \nFruit diameter")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)


c<-ggplot()+
  geom_jitter(data=d.petal,aes(score,pental_lengh_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linespetal,aes(intercept=b_Intercept,slope=b_score),alpha=0.01,color="navy")+
  geom_abline(data=linespetal,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="navy",size=2)+
  geom_abline(data=linespetal.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  geom_abline(data=linespetal.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean \nPetal Length")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)

e<-ggpubr::ggarrange(b,c)

ggpubr::ggarrange(a,e,nrow=2,ncol=1)

save.image("hydro_demand.Rda")

