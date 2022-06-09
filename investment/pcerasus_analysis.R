####FINAL PRUNUS ANAYSIS: see plummywork for scratch
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
#load("plummy.Rda")
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

###Part 1: Turns out phylogeny might matter, or not when we use SE instead of SD

d.sig<-d.flo %>% group_by(specificEpithet) %>% summarise(meanFLS=mean(logFLS),sdFLS=sd(logFLS),nFLS=n(),seFLS=sdFLS / sqrt(nFLS))

###line everthing up for phylosig###
final.df<-d.sig[match(mytree.names, d.sig$specificEpithet),]
namelist2<-final.df$specificEpithet
namelist2==mytree.names
final.df$specificEpithet== mytree.names
#write.csv(final.df,"..//Input/input_clean/FLSdescriptive.csv")

#####

phylosig(pruned.by,final.df$meanFLS,se=final.df$seFLS,method="lambda",nsim = 1000, test=TRUE) #lambda 7.47299e-05 
phylosig(pruned.by,final.df$meanFLS,se=final.df$seFLS,method="K",nsim = 1000, test=TRUE) #K 0.23 


####ordinal model is most descriptive of actual data, so we are going with it here
mod.ord.scale.phlyo<-brm(bbch.v.scale~doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99)) ##

##predict the ordinal
new.data<-data.frame(quant=rep(c( "0%" , "25%",  "50%",  "75%" ,"100%"),13),d.flo%>% group_by(specificEpithet)%>% summarise(doy=quantile(doy)))
new.data$species<-new.data$specificEpithet

season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))



predy2<-fitted(mod.ord.scale.phlyo,newdata = new.data,probs = c(.025,.25,.75,.975))
predy2<-cbind(new.data,predy2)

predy2<-filter(predy2,quant!="100%")
predy2$FLS.scale<-exp(predy2$Estimate)
predy2$Q.25<-exp(predy2$Q25)
predy2$Q.75<-exp(predy2$Q75)

predy3<-predy2 %>%tidyr::gather("phase","likelihood",4:28)
predy3$species2<-predy3$specificEpithet

predy.est<-filter(predy3,str_detect(phase, "^Estimate"))
predy.error<-filter(predy3,str_detect(phase, "^Q2.5"))
predy.error2<-filter(predy3,str_detect(phase, "^Q97.5"))

errorlow <- predy.error %>% 
  group_by(species2,quant) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species2)

errorhigh <- predy.error2 %>% 
  group_by(species2,quant) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species2)

result <- predy.est %>% 
  group_by(species2,quant) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species2)

colnames(errorlow)[3]<-"Q2.5"
colnames(errorhigh)[3]<-"Q97.5"

result$Q2.5<-errorlow$Q2.5
result$Q97.5<-errorhigh$Q97.5
result1<-result

errorlow <- predy.error %>% 
  group_by(species2,quant,phase) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species2)

errorhigh <- predy.error2 %>% 
  group_by(species2,quant,phase) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species2)

result <- predy.est %>% 
  group_by(species2,quant,phase) %>%
  filter(likelihood == max(likelihood)) %>%
  arrange(species2)

colnames(errorlow)[5]<-"Q2.5"
colnames(errorhigh)[5]<-"Q97.5"


result$Q2.5<-errorlow$Q2.5
result$Q97.5<-errorhigh$Q97.5
result
result1$goo<-paste(result1$specificEpithet,result1$phase)
result$goo<-paste(result$specificEpithet,result$phase)

result$most<-NA
result$most<-ifelse(result$goo %in% c(result1$goo),"Y","N")

result$bbch<-NA
result$bbch[which(result$phase=="Estimate.P(Y = 1)" )]<- "BBCH 0"
result$bbch[which(result$phase=="Estimate.P(Y = 2)" )]<- "BBCH 09"
result$bbch[which(result$phase=="Estimate.P(Y = 3)")]<- "BBCH 11"
result$bbch[which(result$phase=="Estimate.P(Y = 4)" )]<- "BBCH 15"
result$bbch[which(result$phase=="Estimate.P(Y = 5)" )]<- "BBCH 17"
result$bbch[which(result$phase=="Estimate.P(Y = 6)" )]<- "BBCH 19"

result$int<-NA
result$int[which(result$phase=="Estimate.P(Y = 1)" )]<- 1
result$int[which(result$phase=="Estimate.P(Y = 2)" )]<- 2
result$int[which(result$phase=="Estimate.P(Y = 3)")]<- 3
result$int[which(result$phase=="Estimate.P(Y = 4)" )]<- 4
result$int[which(result$phase=="Estimate.P(Y = 5)" )]<- 5
result$int[which(result$phase=="Estimate.P(Y = 6)" )]<- 6

result<-filter(result,quant!="100%")
#result<-filter(result,quant!="0%")

### quatify FLS########
makeit<-filter(result,bbch%in% c("BBCH 0","BBCH 09"))
makeit<-makeit[, c('specificEpithet','quant','bbch','likelihood')]

makeit<-tidyr::spread(makeit,bbch,likelihood)

###neeed to do a few alternative descriptors of this for suppliment
makeit$probs<-makeit$`BBCH 0`#+makeit$`BBCH 09`

makeit$cat<-ifelse(makeit$probs>=0.25,"hyserantous","seranthous")
#makeit$cat2<-ifelse(makeit$probs>=0.5,1,0)
#makeit2<-makeit[, c('specificEpithet','quant','cat','cat2')]
result2<-left_join(result,makeit)
result2$likelihood2<-as.numeric(result2$likelihood)

season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))
jpeg("..//Plots/ord_quants_nophylo.jpeg", width=11, height=5,unit="in",res=200)
ggplot(data=result2,aes(bbch,likelihood2))+geom_point()+geom_ribbon(aes(x=int,ymin=0,ymax=likelihood2,fill=cat),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=as.numeric(Q2.5),ymax=as.numeric(Q97.5),width=0))+
  ggthemes::theme_clean(base_size = 10)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  scale_fill_viridis_d()
dev.off()


examplesp<-filter(result2,species2 %in% c("americana","gracilis","maritima","mexicana","subcordata"))

jpeg("..//Plots/ord_quants_exmpsps.jpeg", width=11, height=5,unit="in",res=200)
ggplot(data=examplesp,aes(bbch,likelihood2))+geom_point()+geom_ribbon(aes(x=int,ymin=0,ymax=likelihood2,fill=cat),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=as.numeric(Q2.5),ymax=as.numeric(Q97.5),width=0))+
  ggthemes::theme_clean(base_size = 10)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  scale_fill_viridis_d()
dev.off()

save.image("pcerasus.Rda")
