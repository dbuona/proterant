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

require(mapdata); require(maptools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/Documents/git/proterant/investment/Input")
load("pcerasus.Rda")
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

if(FALSE){
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
pic(final.df$meanFLS,pruned.by,var.contrasts = TRUE,rescaled.tree = TRUE)###
}
####ordinal model is most descriptive of actual data, so we are going with it here
mod.ord.scale.phlyo<-brm(bbch.v.scale~doy+(doy|species)+(doy|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.95,max_treedepth=20)) ##

##predict the ordinal
new.data<-data.frame(quant=rep(c( "0%" , "25%",  "50%",  "75%" ,"100%"),13),d.flo%>% dplyr::group_by(specificEpithet,species)%>% dplyr::summarise(doy=quantile(doy)))


season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))



predy2<-fitted(mod.ord.scale.phlyo,newdata = new.data,probs = c(.025,.25,.75,.975))
predy2<-cbind(new.data,predy2)



predy3<-predy2 %>%tidyr::gather("phase","likelihood",5:40)
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



colnames(errorlow)[6]<-"Q2.5"
colnames(errorhigh)[6]<-"Q97.5"

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

colnames(errorlow)[6]<-"Q2.5"
colnames(errorhigh)[6]<-"Q97.5"


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

makeit$probs2<-as.numeric(makeit$`BBCH 0`)+ as.numeric(makeit$`BBCH 09`)

makeit$cat<-ifelse(makeit$probs>=0.25,"hysteranthous","seranthous")
makeit$cat2<-ifelse(makeit$probs2>=0.5,"hysteranthous","seranthous")
makeit$cat3<-ifelse(makeit$probs2>=0.4,"hysteranthous","seranthous")

makeit$classificationA<-ifelse(makeit$cat=="hysteranthous",1,0)
makeit$classificationB<-ifelse(makeit$cat2=="hysteranthous",1,0)
makeit$classificationC<-ifelse(makeit$cat3=="hysteranthous",1,0)



result2<-left_join(result,makeit)
result2$likelihood2<-as.numeric(result2$likelihood)

season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))
jpeg("..//Plots/ord_quants_phylo.jpeg", width=11, height=11,unit="in",res=200)
main<-ggplot(data=result2,aes(bbch,likelihood2))+geom_point()+geom_ribbon(aes(x=int,ymin=0,ymax=likelihood2,fill=cat2),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=as.numeric(Q2.5),ymax=as.numeric(Q97.5),width=0))+
  ggthemes::theme_clean(base_size = 10)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  scale_fill_viridis_d()+ylab("likelihood")+xlab("vegetative BBCH stage while flowering")+theme(legend.title=element_blank())
dev.off()

alt1<-ggplot(data=result2,aes(bbch,likelihood2))+geom_point()+geom_ribbon(aes(x=int,ymin=0,ymax=likelihood2,fill=cat),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=as.numeric(Q2.5),ymax=as.numeric(Q97.5),width=0))+
  ggthemes::theme_clean(base_size = 10)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  scale_fill_viridis_d()+ylab("likelihood")+xlab("vegetative BBCH stage while flowering")+theme(legend.title=element_blank())

alt2<-ggplot(data=result2,aes(bbch,likelihood2))+geom_point()+geom_ribbon(aes(x=int,ymin=0,ymax=likelihood2,fill=cat3),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=as.numeric(Q2.5),ymax=as.numeric(Q97.5),width=0))+
  ggthemes::theme_clean(base_size = 10)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  scale_fill_viridis_d()+ylab("likelihood")+xlab("vegetative BBCH stage while flowering")+theme(legend.title=element_blank())


jpeg("..//Plots/ord_quants_phylo.jpeg", width=11, height=13,unit="in",res=250)
ggpubr::ggarrange(main,alt1,alt2,ncol=1,nrow=3,common.legend = TRUE,labels=c("a)","b)","c"))
dev.off()

examplesp<-filter(result2,species2 %in% c("americana","angustifolia","maritima","mexicana","subcordata"))

jpeg("..//Plots/ord_quants_exmpsps.jpeg", width=11, height=8,unit="in",res=200)
ggplot(data=examplesp,aes(bbch,likelihood2))+geom_point()+geom_ribbon(aes(x=int,ymin=0,ymax=likelihood2,fill=cat2),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=as.numeric(Q2.5),ymax=as.numeric(Q97.5),width=0))+
  ggthemes::theme_clean(base_size = 10)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  scale_fill_viridis_d()+ylab("likelihood")+xlab("vegetative BBCH stage while flowering")+theme(legend.title=element_blank())
dev.off()


####now model associations
FLSindexA<-dplyr::select(makeit,specificEpithet,quant,classificationA)
FLSindexB<-dplyr::select(makeit,specificEpithet,quant,classificationB)
FLSindexC<-dplyr::select(makeit,specificEpithet,quant,classificationC)

FLSindexA<-tidyr::spread(FLSindexA,quant,classificationA)
FLSindexB<-tidyr::spread(FLSindexB,quant,classificationB)
FLSindexC<-tidyr::spread(FLSindexC,quant,classificationC)

FLSindexA$hystscoreA<-(FLSindexA$`0%`+FLSindexA$`25%`+FLSindexA$`50%`+FLSindexA$`75%`)
FLSindexB$hystscoreB<-(FLSindexB$`0%`+FLSindexB$`25%`+FLSindexB$`50%`+FLSindexB$`75%`)
FLSindexC$hystscoreC<-(FLSindexC$`0%`+FLSindexC$`25%`+FLSindexC$`50%`+FLSindexC$`75%`)


FLSindexA<-dplyr::select(FLSindexA,specificEpithet,hystscoreA)
FLSindexB<-dplyr::select(FLSindexB,specificEpithet,hystscoreB)
FLSindexC<-dplyr::select(FLSindexC,specificEpithet,hystscoreC)

FLSindex<-left_join(FLSindexA,FLSindexB)
FLSindex<-left_join(FLSindex,FLSindexC)


####pdsi model
d<-read.csv("input_clean/pdsi_spi.csv")
d<-dplyr::select(d,-X)
d<-left_join(d,FLSindex)
d$species<-d$specificEpithet

mod.pdsi.phylo<-brm(pdsi~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
mod.pdsi.phyloB<-brm(pdsi~hystscoreB+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
mod.pdsi.phyloC<-brm(pdsi~hystscoreC+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 

chico<-d %>% group_by(specificEpithet) %>% summarise(mean(pdsi,na.rm=TRUE))

exreme<-filter(d, hystscoreB %in% c(0,4))



ggplot(exreme,aes(as.factor(hystscoreB),pdsi))+geom_boxplot()
ggplot(exreme,aes(as.factor(hystscoreB),))+geom_boxplot()
geom_point(aes(color=specificEpithet))
range(d$pdsi,na.rm=TRUE)
mod.min.pdsi.phylo<-brm(pdsi.min~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
##for pdsi, remove phylo since its a species trait not an enviromental trail
#mod.pdsi.nophylo<-brm(pdsi~hystscoreA+(1|specificEpithet),data=d,family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
#mod.pdsi.nophyloB<-brm(pdsi~hystscoreB+(1|specificEpithet),data=d,family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
#mod.pdsi.nophyloC<-brm(pdsi~hystscoreC+(1|specificEpithet),data=d,family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs


#mod.pdsi.nopool<-brm(pdsi~hystscoreA,data=d,family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs

#summary(mod.pdsi.nophylo)
#summary(mod.pdsi.nophyloB)

#summary(mod.pdsi.nophyloC)
#summary(mod.pdsi.nopool)


fixef(mod.pdsi.phyloC,prob=c(.025,.25,.75,.975))[2,]

tab<-data.frame(t(round(fixef(mod.pdsi.phyloB,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tab2<-data.frame(t(round(fixef(mod.pdsi.phylo,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tab3<-data.frame(t(round(fixef(mod.pdsi.phyloC,prob=c(.025,.25,.75,.975))[2,],digits=3)))

tab<-rbind(tab,tab2,tab3)
tab$classification<-c("main analaysis","alternate 1","alternate 2")
tab$Hystanthous_if<-c("50% fl. likelihood  with BBCH 0 & 09","25% fl. likelihood with BBCH 0","40% fl. likelihood with BBCH 0 & 09")
tab$mod_variable<-"mean pdsi"


tabfl<-data.frame(t(round(fixef(mod.petal.phyloB,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfl2<-data.frame(t(round(fixef(mod.petal.phylo,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfl3<-data.frame(t(round(fixef(mod.petal.phyloC,prob=c(.025,.25,.75,.975))[2,],digits=3)))

tabfl<-rbind(tabfl,tabfl2,tabfl3)
tabfl$classification<-c("main analaysis","alternate 1","alternate 2")
tabfl$Hystanthous_if<-c("50% fl. likelihood  with BBCH 0 & 09","25% fl. likelihood with BBCH 0","40% fl. likelihood with BBCH 0 & 09")
tabfl$mod_variable<-"petal length"

tabfr<-data.frame(t(round(fixef(mod.fruit.phyloB,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfr2<-data.frame(t(round(fixef(mod.fruit.phylo,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfr3<-data.frame(t(round(fixef(mod.fruit.phyloC,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfr<-rbind(tabfr,tabfr2,tabfr3)
tabfr$classification<-c("main analaysis","alternate 1","alternate 2")
tabfr$Hystanthous_if<-c("50% fl. likelihood  with BBCH 0 & 09","25% fl. likelihood with BBCH 0","40% fl. likelihood with BBCH 0 & 09")
tabfr$mod_variable<-"fruit diameter"


suptab<-rbind(tab,tabfl,tabfr)
colnames(suptab)
suptab<-suptab[, c(9, 7, 8, 1,2,3,4,5,6)]
xtable::xtable(suptab)


###other covariates
d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")

d.fruit<-filter(d.fruit,fruit_type=="fleshy")

d.petal<-dplyr::select(d.petal,-X)
d.fruit<-dplyr::select(d.fruit,-X)

d.petal<-left_join(d.petal,FLSindex)
d.fruit<-left_join(d.fruit,FLSindex)

d.petal$species<-d.petal$specificEpithet
d.fruit$species<-d.fruit$specificEpithet


mod.petal.phylo<-brm(pental_lengh_mm~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 
mod.petal.phyloB<-brm(pental_lengh_mm~hystscoreB+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99))
mod.petal.phyloC<-brm(pental_lengh_mm~hystscoreC+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99))

fixef(mod.petal.phylo,probs = c(.25,.75))
fixef(mod.petal.phyloB,probs = c(.25,.75))
fixef(mod.petal.phyloC,probs = c(.25,.75))

mod.fruit.phylo<-brm(fruit_diam_mm~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 
mod.fruit.phyloB<-brm(fruit_diam_mm~hystscoreB+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 
mod.fruit.phyloC<-brm(fruit_diam_mm~hystscoreC+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 

fixef(mod.fruit.phylo,probs = c(.25,.75))
fixef(mod.fruit.phyloB,probs = c(.25,.75))
fixef(mod.fruit.phyloC,probs = c(.25,.75))
##B


###plot all that We're choosing scenario B as the best measure of hysteranthy
lines.nophylo<-mod.pdsi.phyloB%>%
  spread_draws(b_Intercept,  b_hystscoreB )

a<-ggplot()+
  geom_jitter(data=d,aes(hystscoreB,pdsi),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_hystscoreB),alpha=0.01,color="skyblue3")+
  geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_hystscoreB)),color="navy",size=2)+
  #geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean \nPDSI at \n collection sites")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)


linespetal<-mod.petal.phyloB%>%
  spread_draws(b_Intercept,  b_hystscoreB )

linesfruit<-mod.fruit.phyloB%>%
  spread_draws(b_Intercept,  b_hystscoreB )
b<-ggplot()+
  geom_jitter(data=d.fruit,aes(hystscoreB,fruit_diam_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linesfruit,aes(intercept=b_Intercept,slope=b_hystscoreB),alpha=0.01,color="skyblue3")+
  geom_abline(data=linesfruit,aes(intercept=mean(b_Intercept),slope=mean(b_hystscoreB)),color="navy",size=2)+
  #geom_abline(data=linesfruit.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=linesfruit.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Fruit diameter")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly \nseason","Through \nmid \nseason","Through \nlate \nseason"))+ggthemes::theme_few(base_size = 11)


c<-ggplot()+
  geom_jitter(data=d.petal,aes(hystscoreB,pental_lengh_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.48,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linespetal,aes(intercept=b_Intercept,slope=b_hystscoreB),alpha=0.01,color="skyblue3")+
  geom_abline(data=linespetal,aes(intercept=mean(b_Intercept),slope=mean(b_hystscoreB)),color="navy",size=2)+
  # geom_abline(data=linespetal.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #  geom_abline(data=linespetal.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Petal Length")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly \nseason","Through \nmid \nseason","Through \nlate \nseason"))+ggthemes::theme_few(base_size = 11)

e<-ggpubr::ggarrange(b,c,labels=c("b)","c)"))


jpeg("..//Plots/dataplots.jpeg", width=8, height=6,unit="in",res=200)
ggpubr::ggarrange(a,e,nrow=2,ncol=1,labels =c("a)" ))
dev.off()




###do species respond plastically?

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

joiner<-dplyr::select(d.um,specificEpithet,lat,year,lon,bbch.v.scale,doy)

pdsi.dater2<-dplyr::left_join(joiner,pdsi.dater)
pdsi.dater2<-pdsi.dater2 %>% distinct()

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

pdsi.dater2$pdsi.z<-zscore(pdsi.dater2$pdsi)
pdsi.dater2$doy.z<-zscore(pdsi.dater2$doy)
pdsi.dater2$species<-pdsi.dater2$specificEpithet
table(pdsi.dater2$species)
pdsi.dater2 %>%group_by(species) %>%summarize(mean(pdsi,na.rm=TRUE),sd(pdsi,na.rm=TRUE))

plastic.mod<-brm(bbch.v.scale~doy.z+pdsi.z+(pdsi.z|specificEpithet)+(1|gr(species, cov = A)),data=pdsi.dater2,data2=list(A=A),family=cumulative("logit"),control = list(adapt_delta=.99),warmup=3000,iter=4000)
plastic.mod.noslp<-brm(bbch.v.scale~doy.z+pdsi.z+(1|specificEpithet)+(1|gr(species, cov = A)),data=pdsi.dater2,data2=list(A=A),family=cumulative("logit"),control = list(adapt_delta=.99),warmup=3000,iter=4000)

fixef(plastic.mod,probs = c(.25,.75))
fixef(plastic.mod.noslp,probs = c(.25,.75))

save.image("pcerasus.Rda")
