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

d.flo$id<-d.flo$specificEpithet ## whoops over wrote the id column here but we dont need if
d.flo$logFLS<-log(d.flo$bbch.v.scale) ## make FLS linear

###Part 1: Turns out phylogeny might matter, or not when we use SE instead of SD

d.sig<-d.flo %>% group_by(specificEpithet) %>% summarise(meanFLS=mean(logFLS),sdFLS=sd(logFLS),nFLS=n(),seFLS=sdFLS / sqrt(nFLS))

###line everthing up for phylosig###
final.df<-d.sig[match(mytree.names, d.sig$specificEpithet),]
namelist2<-final.df$specificEpithet
namelist2==mytree.names
final.df$specificEpithet== mytree.names
#####

phylosig(pruned.by,final.df$meanFLS,se=final.df$seFLS,method="lambda",nsim = 1000, test=TRUE) #lambda 7.47299e-05 
phylosig(pruned.by,final.df$meanFLS,se=final.df$seFLS,method="K",nsim = 1000, test=TRUE) #K 0.23 

####just to be safe we run a phylo mixed modesl
#mod.gaus.scale<-brm(logFLS~doy+(doy|id),data=d.flo,data=d.flo,family=gaussian, control= list(adapt_delta=0.95),warmup = 3000,iter=4000)

#mod.gaus.scale.phylo<-brm(logFLS~doy+(doy|id)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=gaussian, control= list(adapt_delta=0.95),warmup = 3000,iter=4000)

mod.ord.scale.phlyo<-brm(bbch.v.scale~doy+(doy|id)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99)) ## 1 divergent transition

#### now run models for each variable
d.pdsi<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi.csv")
#d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
#d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")
#d.phen<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruit_phen.csv")
#d.cold<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi_wint.csv")
#d.fruit<-filter(d.fruit,fruit_type=="fleshy")


#z.score everything for future analyses
zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
d.pdsi$pdsi.z<-zscore(d.pdsi$pdsi)
#d.petal$petal.z<-zscore(d.petal$pental_lengh_mm)
#d.fruit$fruit.z<-zscore(d.fruit$fruit_diam_mm)
#d.phen$phen.z<-zscore(d.phen$doy)


### get ready to play with phylogeny
d.pdsi<-dplyr::filter(d.pdsi,specificEpithet %in% pruned.by$tip.label)
d.pdsi$ID<-d.pdsi$specificEpithet

#d.petal$ID<-d.petal$specificEpithet
#d.fruit$ID<-d.fruit$specificEpithet
pdsi.mod.phylo<-brm(pdsi~(1|ID)+(1|gr(specificEpithet, cov = A)),data=d.pdsi,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.98))

##add min pdsi
#d.pdsi$pdsi.min.z<-zscore(d.pdsi$pdsi.min)

if (FALSE){ ## if you decide there is too much colinearity you can use these
  pdsi.mod.z<-brm(pdsi.z~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
  pdsi.min.mod.z<-brm(pdsi.min.z~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
  petal.mod.z<- brm(petal.z~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
  fruit.mod.z<- brm(fruit.z~(1|id)+(1|specificEpithet),data=d.fruit,warmup=3000,iter=4000)
pdsi.mod.phylo<-brm(pdsi~(1|ID)+(1|gr(specificEpithet, cov = A)),data=d.pdsi,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.98))
petalmod.phylo<- brm(pental_lengh_mm~(1|ID)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.petal,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.98))
fruitmod.phylo<- brm(fruit_diam_mm~(1|ID)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.fruit,data2=list(A=A),warmup=4500,iter=6000,control=list(adapt_delta=0.98)) # 4 divergent transitions

}


pdsi.mod<-brm(pdsi~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
petal.mod<- brm(pental_lengh_mm~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
fruit.mod<- brm(fruit_diam_mm~(1|id)+(1|specificEpithet),data=d.fruit,warmup=2500,iter=4000)

if (FALSE){
pdsi.mod.z.phylo<-brm(pdsi.z~(1|ID)+(1|gr(specificEpithet, cov = A)),data=d.pdsi,data2=list(A=A),warmup=4000,iter=5000,control=list(adapt_delta=0.98))
petalmod.z.phylo<- brm(petal.z~(1|ID)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.petal,data2=list(A=A),warmup=4000,iter=5000,control=list(adapt_delta=0.99)) # 1 dierget transition
fruitmod.z.phylo<- brm(fruit.z~(1|ID)+(1|id)+(1|gr(specificEpithet, cov = A)),data=d.fruit,data2=list(A=A),warmup=4000,iter=5000,control=list(adapt_delta=0.99))
}
summary(petalmod.z.phylo)


pdsiout<-dplyr::select(as.data.frame(coef(pdsi.mod.z)),1:2)
pdsiminout<-dplyr::select(as.data.frame(coef(pdsi.min.mod.z)),1:2)
petalout<-dplyr::select(as.data.frame(coef(petal.mod.z)$specificEpithet),1:2)
fruitout<-dplyr::select(as.data.frame(coef(fruit.mod.z)$specificEpithet),1:2)
FLSout<-dplyr::select(as.data.frame(coef(mod.gaus.scale.phylo)$specificEpithet),1:2)
#phenout<-dplyr::select(as.data.frame(coef(phenmod.z)),1:2)

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

cor(newdat$pdsi_mean,newdat$petal_mean, use="pairwise.complete.obs") #0.542
cor(newdat$fruit_mean,newdat$petal_mean, use="pairwise.complete.obs") #0.493
cor(newdat$fruit_mean,newdat$pdsi_mean, use="pairwise.complete.obs") # 0.547



save.image("plummy.Rda")


##predict the ordinal
new.data<-data.frame(quant=rep(c( "0%" , "25%",  "50%",  "75%" ,"100%"),13),d.flo%>% group_by(specificEpithet)%>% summarise(doy=quantile(doy)))
new.data$id<-new.data$specificEpithet


predy<-fitted(mod.gaus.scale.phylo,newdata = new.data,probs = c(.25,.75))
predy<-cbind(new.data,predy)

season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))

predy<-filter(predy,quant!="100%")
predy$FLS.scale<-exp(predy$Estimate)
predy$Q.25<-exp(predy$Q25)
predy$Q.75<-exp(predy$Q75)
season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))
ggplot(data=predy,aes(quant,FLS.scale))+geom_point(aes(color=id))+
  geom_errorbar(aes(ymin=Q.25,ymax=Q.75,width=0,color=id))+ggthemes::theme_clean(base_size = 11)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+
  theme(strip.text = element_text(face = "italic"))+
  scale_y_continuous(breaks=c(1,2,3,4,5),labels = c("BBCH 0","BBCH 09","BBCH 11","BBCH 15","BBCH 17"),name = "BBCH")

dev.off()

predy2<-fitted(mod.ord.scale.phlyo,newdata = new.data,probs = c(.025,.975))
predy2<-cbind(new.data,predy2)



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
makeit$probs<-makeit$`BBCH 0`#+makeit$`BBCH 09`

makeit$cat<-ifelse(makeit$probs>=0.25,"hyserantous","seranthous")
#makeit$cat2<-ifelse(makeit$probs>=0.5,1,0)
#makeit2<-makeit[, c('specificEpithet','quant','cat','cat2')]
result2<-left_join(result,makeit)
result2$likelihood2<-as.numeric(result2$likelihood)

season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))
jpeg("..//Plots/ord_quants.jpeg", width=11, height=5,unit="in",res=200)
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
#always1
sos<-c(1,1,1,1,1,1,1,1,1,0,0,0,1)
es<-c(0,0,1,1,0,1,1,0,0,0,0,0,1)
ms<-c(0,0,1,0,0,1,1,0,0,0,0,0,1)
ls<-c(0,0,0,0,0,0,1,0,0,0,0,0,1)


hystscore<-data.frame(specificEpithet=sort(unique(d.flo$specificEpithet)),sos,es,ms,ls)
hystscore$score <- rowSums(hystscore[2:5])
write.csv(hystscore,"input_clean/hystscore.csv")

### make tree
newdat<-left_join(newdat,hystscore)

tree2<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre")


names.intree<-tree2$tip.label # names the names
namelist<-unique(d.flo$specificEpithet)

to.prune<-which(!names.intree%in%namelist) #prun the tree
pruned.by2<-drop.tip(tree,to.prune)
plotTree(pruned.by2)# this is the tree


###what are the tip labels in pruned phylogeny?

mytree.names<-pruned.by2$tip.label # did i get them all
intersect(namelist,mytree.names) #yes







### Q1 DO we bother with phylogeny? Only 7 species, and as you can see, weak signal
meanfls<-d.flo %>% dplyr::group_by(specificEpithet)%>% dplyr::summarise(meanFLS=mean(bbch.v.scale),sdFLS=sd(bbch.v.scale))##Take mean FLS values
d.phylo<-filter(d.flo,specificEpithet %in% mytree.names)

d.phylo<-left_join(d.phylo,meanfls)
##quick phylo sig
meanfls4phylo<-filter(meanfls,specificEpithet %in% mytree.names)
meanfls4phylo<-left_join(meanfls4phylo,hystscore)




library(ggtree)

meanfls4phylo$score2<-as.factor(meanfls4phylo$score)
jpeg("..//Plots/phylosig2.jpeg", width=4, height=4,unit="in",res=300)
p<-ggtree(pruned.by2,layout = "roundrect")
p %<+% meanfls4phylo+geom_tiplab(hjust=-.2,align=TRUE,fontface="italic")+geom_tippoint(aes(color=score2),size=5,shape=15)+ xlim(0, 3)+scale_color_viridis_d(option="turbo",name="Hysteranthy",  labels=c("Never","At start of season","Through early season","Through mid season","Through late season"))
dev.off()

##full model
head(newdat)
newdat$interscore<-newdat$score+1
fullymody<-brm(interscore~me(pdsi_mean,pdsi_se)+me(petal_mean,petal_se)+me(fruit_mean,fruit_se),data=newdat,warmup=3000,iter=4000, family=cumulative("logit"))
fixef(fullymody)
####model each
catgaus.pdsi<-brm(pdsi_mean|mi(pdsi_se)~score,data=newdat,family=gaussian(),control = list(adapt_delta=.99),warmup=5000,iter=6000)
fixef(catgaus.pdsi)


fixef(catgaus.pdsi)
get_variables(catgaus.pdsi)
lines<-catgaus.pdsi%>%
  spread_draws(b_Intercept,  b_score )

catgaus.fruit<-brm(fruit_mean|mi(fruit_se)~score,data=newdat,family=gaussian(),control = list(adapt_delta=.99),warmup=5000,iter=6000)
fixef(catgaus.fruit)
catgaus.fruit<-add_criterion(catgaus.fruit,"loo")

lines2<-catgaus.fruit%>%
  spread_draws(b_Intercept,  b_score )

catgaus.petal<-brm(petal_mean|mi(petal_se)~score,data=newdat,family=gaussian(),control = list(adapt_delta=.99),warmup=5000,iter=6000)
fixef(catgaus.petal)
catgaus.petal<-add_criterion(catgaus.petal,"loo")

loo_compare(catgaus.petal,catgaus.fruit)

lines3<-catgaus.petal%>%
  spread_draws(b_Intercept,  b_score )

##extrac and group posteriors I think these are doing raneffs
goober2<-pdsi.mod.phylo%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(goober2)

flooby<-left_join(goober2,hystscore)
flooby<-flooby%>% group_by(score)
flooby$Intercept<-flooby$Intercept+0.16936260


###peta
gooberp<-petal.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 

flooby2<-left_join(gooberp,hystscore)
flooby2<-flooby2%>% group_by(score)
flooby2$Intercept<-flooby2$Intercept+5.0609407

###fruit
gooberf<-fruit.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 

flooby3<-left_join(gooberf,hystscore)
flooby3<-flooby3%>% group_by(score)
flooby3$Intercept<-flooby3$Intercept+19.22416



gooberx<-pdsi.min.mod.z%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(gooberx)

floobyx<-left_join(gooberx,hystscore)
floobyx<-floobyx%>% group_by(score)
floobyx$Intercept<-floobyx$Intercept+0.16936260


ggplot()+
  #geom_abline(data=lines,aes(intercept=b_Intercept,slope=b_score),alpha=0.009,color="navy")+
  #geom_abline(data=lines,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="black")+
  stat_eye(data=floobyx,aes(score,Intercept),alpha=0.6,fill="grey50")+
  ylab("Mean \nPDSI at \n collection sites")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)


a<-ggplot()+
  #geom_abline(data=lines,aes(intercept=b_Intercept,slope=b_score),alpha=0.009,color="navy")+
  #geom_abline(data=lines,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="black")+
  stat_eye(data=flooby,aes(score,Intercept),alpha=0.6,fill="grey50")+
 ylab("Mean \nPDSI at \n collection sites")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                   labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)



b<-ggplot()+
  geom_abline(data=lines2,aes(intercept=b_Intercept,slope=b_score),alpha=0.009,color="navy")+
  geom_abline(data=lines2,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="black")+
  stat_eye(data=flooby3,aes(score,Intercept),alpha=0.6,fill="grey50")+

ylab("Mean \nfruit diameter")+xlab("FLS group")+scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                                                                                                         labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)
c<-ggplot()+
  geom_abline(data=lines3,aes(intercept=b_Intercept,slope=b_score),alpha=0.009,color="navy")+
geom_abline(data=lines3,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="black")+
stat_eye(data=flooby2,aes(score,Intercept),alpha=0.6,fill="grey50")+
ylab("Mean \npetal length")+xlab("FLS group")+scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                                                                                                                    labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)




d<-ggpubr::ggarrange(b,c,nrow=1,ncol=2,labels=c("b)","c)"))
e<-ggpubr::ggarrange(a,d,nrow=2,ncol=1,labels=c("a)",""))

jpeg("..//Plots/dataplots.jpeg", width=7, height=5,unit="in",res=300)
e
dev.off()

###plastic
d.um<-d.flo
palmer.b <- brick("..//Data/lbda-v2_kddm_pmdi_2017.nc")

lonpoints<-d.um$lon # make vector of prunus coordinates
latpoints<-d.um$lat #
extract.pts <- cbind(lonpoints,latpoints)
palmer.b
palmer.b2 <-palmer.b[[1800:2018]]## subset to pnly last century
palmer.b2<-brick(palmer.b2)
ext<-raster::extract(palmer.b2,extract.pts,method="simple")
colnames(ext)<-(c(1799:2017))
ext<-as.data.frame(ext)
ext$lat<-latpoints
ext$lon<-lonpoints
colnames(ext)
pdsi.dater<-tidyr::gather(ext,"year","pdsi",1:219)
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

moda<-brm(bbch.v.scale~pdsi.z+doy.z+(pdsi.z+doy.z|specificEpithet),data=pdsi.dater2,family=cumulative("logit"),control = list(adapt_delta=.98))
summary(moda)
+(1|gr(specificEpithet, cov = A)),data=d.um, data2=list(A = A),family=cumulative("logit"))
(doy|ID)+(1|gr(specificEpithet, cov = A)),data=d.flo,
fixef(moda,probs = c(.25,.75))
coef(moda,probs = c(.25,.75))
goo<-as.data.frame(coef(moda))

goo<-dplyr::select(goo,1:4)

goo2<-as.data.frame(fixef(moda,probs =c(.25,.75)))
goo2$species<-rownames(goo2)
goo2<-filter(goo2,species=="pdsi")
goo$species<-rownames(goo)
colnames(goo)<-colnames(goo2)



goo<-rbind(goo,goo2)
goo$species<-ifelse(goo$species=="pdsi","Main Effect",goo$species)
goo$species
goo$Estimate2<-(-goo$Estimate)
goo$Q252<-(-goo$Q25)
goo$Q752<-(-goo$Q75)


get_variables(moda)
coef(moda)
yaya<-moda%>%
  spread_draws(r_specificEpithet[specificEpithet,condition])
yaya<-filter(yaya,condition=="pdsi.z")

yaya2<-moda%>%
  spread_draws(b_pdsi.z)
yaya2$specificEpithet<-"Main effect"

jpeg("..//Plots/droughtstuff.jpg", width=11, height=9,unit="in",res=300)
ggplot()+
  stat_eye(data=yaya,aes(r_specificEpithet,specificEpithet),alpha=.3)
+
  stat_eye(data=yaya2,aes(b_pdsi.z,specificEpithet))
  geom_vline(xintercept=0,linetype="dashed")+scale_y_discrete(name ="species", 
                                                              limits=c("alleghaniensis", "americana"    ,  "angustifolia"  , "gracilis"   ,    "hortulana"    , 
                                                                       "maritima"  ,     "mexicana"    ,   "munsoniana"   ,  "nigra"     ,     "rivularis"   ,  
                                                                       "subcordata" ,    "texana"     ,    "umbellata", "Main effect")) +ggthemes::theme_few(base_size = 11)+
  theme(axis.text.y = element_text(face=ifelse(goo$species=="Main Effect","bold","italic")))+xlab("Estimate")+
  annotate(geom="text",color="black", x=-1, y=13.5,label="Increased aridity decreases \n likelihood \nof flowering-first")+
  annotate(geom="text",color="black", x=1, y=13.5,label="Increased aridity increases \n likelihood \nof flowering-first")+
  annotate("segment", x = .5, xend = 1.5, y = 14.3, yend = 14.3,
           arrow = arrow(ends = "last" , length = unit(.2,"cm")))+
  annotate("segment", x = -.5, xend = -1.5, y = 14.3, yend = 14.3,
           arrow = arrow(ends = "last" , length = unit(.2,"cm")))
dev.off()
?arrow()

goo$effect<-ifelse(goo$species=="Main Effect","main","species")
q<-ggplot(goo,aes(Estimate2,species))+geom_point(aes(size=effect))+
  geom_errorbarh(aes(xmin=Q252,xmax=Q752,height=0))+scale_size_manual(values=c(4,2))+
    geom_vline(xintercept = 0, color="red")+ggthemes::theme_clean(base_size = 11) + theme(axis.text.y = element_text(face=ifelse(goo$species=="Main Effect","bold","italic")))+
  +theme(legend.position = "none")+xlab("Drought effect estimate")+
  annotate(geom="text",color="gray39", x=-1, y=13.5,label="Increased aridity increases \n likelihood \nof flowering-first")+
  annotate(geom="text",color="gray39", x=1, y=13.5,label="Increased aridity decreases \n likelihood \nof flowering-first")+xlim(-1.5,1.5)

jpeg("..//Plots/droughtstuff.jpg", width=11, height=9,unit="in",res=300)
#ggpubr::ggarrange(a,q,nrow=2,heights = c(.5,.7),labels = c("a)","b)"))
q
dev.off()


###quick model to see if pdsi impact day of flowering in general
modcont<-brm(doy.cent~pdsi+lat+(pdsi|specificEpithet),data=d.um)


fixef(modcont,probs = c(.25,.75))
coef(modcont,probs = c(.25,.75))

jpeg("..//Plots/alt_hypothesis.jpg", width=11, height=6,unit="in",res=300)
ggpubr::ggarrange(b,c,d,labels = c("a)","b)","c)"),ncol=3)
dev.off()
####phylogeny


###make a range map

mid.herb<-read.csv("..//SymbOutput_2020-10-26_133412_DwC-A/occurrences.csv")
pruno.ref<-filter(mid.herb,!is.na(decimalLatitude))
table(pruno.ref$stateProvince)
table(pruno.ref$specificEpithet)
pruno.ref$lon<-pruno.ref$decimalLongitude
pruno.ref$lat<-pruno.ref$decimalLatitude

##now add county level coordiantes
#### two without coordinates
pruno.unref<-filter(mid.herb,is.na(decimalLatitude)) ## filter entries with no lat/lon
### prep county data
geoCounty$rMapState<-str_to_title(geoCounty$rMapState) ### centriod coordinates for all US counties
colnames(geoCounty)[6]<-"stateProvince" 

colnames(geoCounty)[c(2,7)]<-c("County_name","county") ## MATCH NAMES

pruno.unref<-filter(pruno.unref,id!=17453277) #remvoe problematic canadian entry
pruno.unref$county<-tolower(pruno.unref$county) ## make ours lower case



pruno.unref$county<-gsub("county","",pruno.unref$county) ### get rid of extenious coounty info
pruno.unref$county<-gsub("()","",pruno.unref$county) # ""
pruno.unref$county<-gsub("co.","",pruno.unref$county,fixed = TRUE) # ""




pruno.unref<-dplyr::left_join(pruno.unref,geoCounty,by=c("county","stateProvince"))
pruno.unref<-dplyr::select(pruno.unref,-County_name,-fips, -state)
intersect(colnames(pruno.unref),colnames(pruno.ref))

pruno.ref$geoclass<-"geo_referenced"
pruno.unref$geoclass<-"county_centroid"
pruno.unref<-filter(pruno.unref,!is.na(lon)) # select only entries with lat lon
pruneo<-rbind(pruno.ref,pruno.unref)



pruneo<- pruneo %>%filter(specificEpithet %in% unique(d.flo$specificEpithet))
pruneo<-filter(pruneo,lon<(-60))
pruneo<-filter(pruneo,lon>(-150))

hull_cyl <- pruneo %>%
  group_by(specificEpithet) %>%
  slice(chull(lon,lat))
colnames(hull_cyl)

hull_cyl<-filter(hull_cyl, specificEpithet %in% c("subcordata","texana")) 

usa <- map_data("usa")


jpeg("..//Plots/map.jpeg", width=11, height=7,unit="in",res=300)
ggplot()+geom_polygon(data=usa,aes(long,lat,group=group),fill="white",color="black")+
geom_point(data=pruneo,aes(lon,lat,color=specificEpithet,shape=specificEpithet),alpha=0.9)+ggthemes::theme_few()+
scale_color_viridis_d(option="turbo")+
  scale_shape_manual(values = rep(c(0,1,2,15,16,17), len = 14))
dev.off()

