####to make potential prunus manuscript figures

#Note  I think there is something missing in how i extract the posteriors. iving the equivelent of ranef instead of coef

#also thsi code time travels (things from below are used to make things above) so if i ever decide to re-run everything its going to need some tweaks
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
library(lubridate)
library(stringr)
library("tidybayes")
library(raster)

require(mapdata); require(maptools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/Documents/git/proterant/investment/Input")
load("predictory.Rda")



##read in cleaned data
d.flo<-read.csv("input_clean/FLS_clean.csv") ##data
tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre") ##tree

##### give tree branch lengths
is.ultrametric(tree)
tree<-compute.brlen(tree, method = "Grafen")## make ultrametric
#check names
names.intree<-tree$tip.label # names the names
namelist<-unique(d.flo$specificEpithet)
to.prune<-which(!names.intree%in%namelist) #prune the tree
pruned.by<-drop.tip(tree,to.prune)
plotTree(pruned.by)# this is the tree

###what are the tip labels in pruned phylogeny?

mytree.names<-pruned.by$tip.label # did i get them all
intersect(namelist,mytree.names) #yes

A <- ape::vcv.phylo(pruned.by) ## make acovarience matrix for brms models

d.flo$species<-d.flo$specificEpithet ## whoops over wrote the id column here but we dont need if
d.flo$logFLS<-log(d.flo$bbch.v.scale) ## make FLS linear



#### now run models for each variable
d.pdsi<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi.csv")
d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")
#d.phen<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruit_phen.csv")
#d.cold<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi_wint.csv")
d.fruit<-filter(d.fruit,fruit_type=="fleshy")

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
d.pdsi$pdsi.z<-zscore(d.pdsi$pdsi)
d.petal$petal.z<-zscore(d.petal$pental_lengh_mm)
d.fruit$fruit.z<-zscore(d.fruit$fruit_diam_mm)
#d.phen$phen.z<-zscore(d.phen$doy)


### get ready to play with phylogeny
d.pdsi<-dplyr::filter(d.pdsi,specificEpithet %in% pruned.by$tip.label)
d.pdsi$species<-d.pdsi$specificEpithet

d.petal$species<-d.petal$specificEpithet
d.fruit$species<-d.fruit$specificEpithet

#d.pdsi$pdsi.min.z<-zscore(d.pdsi$pdsi.min)

 ## if you decide there is too much colinearity you can use these
  pdsi.mod.z<-brm(pdsi.z~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
  #pdsi.min.mod.z<-brm(pdsi.min.z~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
  petal.mod.z<- brm(petal.z~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
  fruit.mod.z<- brm(fruit.z~(1|id)+(1|specificEpithet),data=d.fruit,warmup=3000,iter=4000) #2 divergent
 
  pdsi.mod.phylo<-brm(pdsi~(1|species)+(1|gr(specificEpithet, cov = A)),data=d.pdsi,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.99))
  petalmod.phylo<- brm(pental_lengh_mm~(1|species)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.petal,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.99))
  fruitmod.phylo<- brm(fruit_diam_mm~(1|species)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.fruit,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.99))  
  

pdsi.mod<-brm(pdsi~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
petal.mod<- brm(pental_lengh_mm~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
fruit.mod<- brm(fruit_diam_mm~(1|id)+(1|specificEpithet),data=d.fruit,warmup=2500,iter=4000)


  pdsi.mod.z.phylo<-brm(pdsi.z~(1|species)+(1|gr(specificEpithet, cov = A)),data=d.pdsi,data2=list(A=A),warmup=4000,iter=5000,control=list(adapt_delta=0.99))
  petalmod.z.phylo<- brm(petal.z~(1|species)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.petal,data2=list(A=A),warmup=4000,iter=5000,control=list(adapt_delta=0.99)) # 
  fruitmod.z.phylo<- brm(fruit.z~(1|species)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.fruit,data2=list(A=A),warmup=4000,iter=5000,control=list(adapt_delta=0.99)) # 2 divergent transition

save.image("predictory.Rda")
summary(petalmod.z.phylo)


pdsiout<-dplyr::select(as.data.frame(coef(pdsi.mod)),1:2)
petalout<-dplyr::select(as.data.frame(coef(petal.mod)$specificEpithet),1:2)
fruitout<-dplyr::select(as.data.frame(coef(fruit.mod)$specificEpithet),1:2)



###plot all this
get_variables(pdsi.mod.z)
pdsi.z.plot<-pdsi.mod.z%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

pdsi.plot<-pdsi.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

pdsi.z.plot.phylo<-pdsi.mod.z.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

pdsi.plot.phylo<-pdsi.mod.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

pdsia<-ggplot(pdsi.plot,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()
pdsib<-ggplot(pdsi.z.plot,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()

pdsic<-ggplot(pdsi.plot.phylo,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()
pdsid<-ggplot(pdsi.z.plot.phylo,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()

ggpubr::ggarrange(pdsia,pdsib,pdsic,pdsid)


flower.z.plot<-petal.mod.z%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

flower.plot<-petal.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

flower.z.plot.phylo<-petalmod.z.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

flower.plot.phylo<-petalmod.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])


floa<-ggplot(flower.plot,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()
flob<-ggplot(flower.z.plot,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()

floc<-ggplot(flower.plot.phylo,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()
flod<-ggplot(flower.z.plot.phylo,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()

ggpubr::ggarrange(floa,flob,floc,flod)




fruit.z.plot<-fruit.mod.z%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

fruit.plot<-fruit.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

fruit.z.plot.phylo<-fruitmod.z.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])

fruit.plot.phylo<-fruitmod.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,Estimate])


frua<-ggplot(fruit.plot,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()
frub<-ggplot(fruit.z.plot,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()

fruc<-ggplot(fruit.plot.phylo,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()
frud<-ggplot(fruit.z.plot.phylo,aes(reorder(specificEpithet, -r_specificEpithet),r_specificEpithet))+stat_eye()


ggpubr::ggarrange(frua,frub,fruc,frud)


colnames(pdsiout)<-c("pdsi_mean","pdsi_se")
colnames(pdsiminout)<-c("pdsimin_mean","pdsimin_se")
colnames(petalout)<-c("petal_mean","petal_se")
colnames(fruitout)<-c("fruit_mean","fruit_se")
colnames(FLSout)<-c("FLS_mean","FLS_se")
#colnames(phenout)<-c("phen_mean","phen_se")

#phenout$specificEpithet<-rownames(phenout)
fruitout$specificEpithet<-rownames(fruitout)
petalout$specificEpithet<-rownames(petalout)
pdsiout$specificEpithet<-rownames(pdsiout)
pdsiminout$specificEpithet<-rownames(pdsiminout)
FLSout$specificEpithet<-rownames(FLSout)

newdat<-left_join(fruitout,pdsiout)
newdat<-left_join(newdat,petalout)
newdat<-left_join(newdat,FLSout)

