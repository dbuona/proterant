##American plums 4KNB
####FINAL PRUNUS ANAYSIS: see plummywork for scratch and hydraulic demand
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
graphics.off()
library(dplyr)
library(ggplot2)
library(rstan)
library(brms)
library(phytools)
library(ape)
library(tidybayes)



setwd("~/Documents/git/proterant/investment/Input")
#load("..//forKNB/amer_plums_KNB.Rda")
d<-read.csv("..//forKNB/amerplum_final.csv") ##read in data
tree<-read.tree("~/Documents/git/proterant/investment/forKNB/amerplum_tree.tre") ##read in tree

A <- ape::vcv.phylo(tree) ### build VCv matrix based on phylogeny

###main model####
mod.ord.<-brm(bbch.v.scale~YEAR.hin+doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
##################

###model w/o doy for Supporting Infortmation####
mod.ord.nodoy<-brm(bbch.v.scale~YEAR.hin+(1|species)+(1|gr(specificEpithet, cov = A)),data=d,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
######################

####################start to make Figure 2##############
#make new data frame for predictions#########
season<-d%>% dplyr::group_by(specificEpithet,species) %>% summarise(start=min(doy),end=max(doy))
#lat<-d.flo%>% dplyr::group_by(specificEpithet,species) %>% summarise(lat=mean(lat,na.rm=TRUE))
sps<-season$species

new.dat<-data.frame()
for (z in c(1:length(sps))){
  dat <-filter(season,species==sps[z])
  goo<-data.frame(doy=seq(dat$start,dat$end,by=1),species=rep(sps[z],length(dat$start)))### this establshed the length of the flowering season for each species
  new.dat<-rbind(goo,new.dat)
}

new.dat$YEAR.hin<-20 ### since Year doesn't effect estimates we choose a standard year  (2000) for all species predictors
new.dat$specificEpithet<-new.dat$species #### add the phylogeny variable
###########################################

#####make posterior predictions
review_pred <- mod.ord. %>% 
  epred_draws(newdata =new.dat,ndraws=1000,seed = 1234)

#### classify hysteranthy as BBCH 0,9,11
review_pred$class<-ifelse(review_pred$.category %in% c(1,2,3),"hysteranthous","non-hysteranthous")

#### classify hysteranthy as BBCH 0,9,
review_pred$class2<-ifelse(review_pred$.category %in% c(1,2),"hysteranthous","non-hysteranthous")

#### group for an average species level prediction of hysteranthy per species for each day of season with both
###to make hysteranthy index
check<-review_pred %>% group_by(species,class,.draw,doy)%>%summarise(pred=sum(.epred))
check2<-review_pred %>% group_by(species,class2,.draw,doy)%>%summarise(pred=sum(.epred))

#### make figre 2a
cofs<-as.data.frame(fixef(mod.ord.,probs = c(.055,.25,.75,.945)))
cofs<-tibble::rownames_to_column(cofs,"predictor")
cofs<-filter(cofs,predictor %in% c("doy","YEAR.hin"))
p2a<-ggplot(cofs,aes(Estimate,predictor))+geom_point(size=4)+
  geom_errorbarh(aes(xmin=`Q5.5`,xmax=`Q94.5`),height=0,size=.5)+
  geom_errorbarh(aes(xmin=`Q25`,xmax=`Q75`),height=0,size=1)+geom_vline(xintercept=0)+
  scale_y_discrete(labels=c("day of year","year"))+xlab("Estimated effect")+ggthemes::theme_few()


#### make index data per species
sumz<-check %>%group_by(species,class)%>% summarise(index=round(mean(pred),2))
sumz<-filter(sumz,class=="hysteranthous")

sumz2<-check2 %>%group_by(species,class2)%>% summarise(index=round(mean(pred),2))
sumz2<-filter(sumz2,class2=="hysteranthous")

sumz<-arrange(sumz,-index)#### this is a list of each species hysteranthy index score
sumz2<-arrange(sumz2,-index) #### this is a list of each species hysteranthy index score with alternate metric

#Make figure 2b
p2b<-ggplot(check,aes(doy,pred))+geom_line(aes(color=class,group=.draw),size=0.3)+facet_wrap(~factor(species,levels=c("mexicana","umbellata","angustifolia","maritima","gracilis","americana","munsoniana","alleghaniensis","nigra","hortulana","texana","rivularis","subcordata")))+
  ggthemes::theme_few(base_size =11)+
  scale_color_viridis_d(direction = -1,begin=.6,end=.1)+geom_label(data=sumz,x=150,y=.75,aes(label=index))+
  xlab("day of year")+
  ylab("hysteranthy likelihood")+
  theme(strip.text = element_text(face = "italic"))+theme(legend.position = "bottom",legend.title = element_blank())

#####make supplimental figures
##S1
ggplot(d,aes(doy))+geom_histogram(bins=15)+facet_wrap(~species,scales="free_x")+
  ggthemes::theme_few()+theme(strip.text = element_text(face="italic"))+xlab("day of year collected")

#S2
###predict from model
simpdat<-data.frame(species=unique(d$specificEpithet),specificEpithet=unique(d$specificEpithet))
simpdat$YEAR.hin<-20
prezsim<-fitted(mod.ord.nodoy,newdata=simpdat,probs = c(.055,.955))
prezsim<-cbind(simpdat,prezsim)

####pull mean estimates
prezest<- prezsim %>%
  dplyr::select(species,contains("Estimate"))
colnames(prezest)[2:7]<-c(1,2,3,4,5,6)
prezest$index.nodoy<-prezest$`1`+prezest$`2`+prezest$`3`
prezest$index2.nodoy<-prezest$`1`+prezest$`2`
prezest$index3.nodoy<-prezest$`1`
prezest<-tidyr::gather(prezest,"stage","probability",2:7)
###each index is a different classifiaction of hysteranthy (bbch 0,9,11, or bbch 0,9 or bbch 0 only)

###low bound CI
prezlow<- prezsim %>%
  dplyr::select(species,contains("Q5"))
colnames(prezlow)[2:7]<-c(1,2,3,4,5,6)
prezlow<-tidyr::gather(prezlow,"stage","low",2:7)

###upper bound CI
prezhigh<- prezsim %>%
  dplyr::select(species,contains("Q9"))
colnames(prezhigh)[2:7]<-c(1,2,3,4,5,6)
prezhigh<-tidyr::gather(prezhigh,"stage","high",2:7)

###merge
prezest<-left_join(prezest,prezlow)
prezest<-left_join(prezest,prezhigh)
prezest<-distinct(prezest)
prezest$stage<-as.numeric(prezest$stage)

#### clasify hysteranthy
prezest$class<-ifelse(prezest$stage <4,"hysteranthous","non-hysteranthous")
indy<-dplyr::select(prezest,species,index.nodoy,index2.nodoy,index3.nodoy)
indy<-distinct(indy)
indy$index.nodoy<-round(indy$index.nodoy,2)

##relabel numeric stages back to orginal BBCH
prezest$bbch<-NA
prezest$bbch[which(prezest$stage==1 )]<- "BBCH 0"
prezest$bbch[which(prezest$stage==2 )]<- "BBCH 09"
prezest$bbch[which(prezest$stage==3)]<- "BBCH 11"
prezest$bbch[which(prezest$stage==4)]<- "BBCH 15"
prezest$bbch[which(prezest$stage==5 )]<- "BBCH 17"
prezest$bbch[which(prezest$stage==6 )]<- "BBCH 19"


#plot
ggplot(prezest,aes(bbch,probability))+
  geom_point(size=1)+geom_ribbon(aes(x=as.numeric(stage),ymin=0,ymax=probability,fill=class,alpha=class))+
  facet_wrap(~factor(species,levels=c("mexicana","umbellata","angustifolia","maritima","gracilis","munsoniana","alleghaniensis","nigra","americana","texana","rivularis","hortulana","subcordata")))+
  geom_errorbar(aes(ymin=low,ymax=high),width=0)+
  theme(strip.text = element_text(face = "italic"))+
  geom_label(data=indy,x=5,y=.5,aes(label=index.nodoy))+
  scale_alpha_manual(values=c(0.6,0))+ggthemes::theme_few()+
  scale_fill_viridis_d(direction = -1,begin=.6,end=.1)+
  theme(strip.text = element_text(face = "italic"))+xlab("")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 270,vjust =.4,size = 7))


##read in average trait values for the species
spsmean<-read.csv("..//forKNB/speciestraits.csv")

###join to hysteranthy data
sumz<-left_join(sumz,spsmean)

##make zscoring function
zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

sumz$pdsi.z<-zscore(sumz$mean.pdsi)
sumz$petal.z<-zscore(sumz$mean.petal)
sumz$temp.z<-zscore(sumz$meanspringT)

##check covarience
cor(sumz$petal.z,sumz$temp.z) ### -0.72 maybe too high
cor(sumz$petal.z,sumz$pdsi.z) ##-0.1 ok

mod.review.PDSI<- brms::brm(
  brms::bf(index ~ pdsi.z*petal.z,
           phi ~1),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 5000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

mod.review.Temp<- brms::brm(
  brms::bf(index ~ temp.z*petal.z,
           phi ~1),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 5000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

bayes_R2(mod.review.PDSI,probs = c(.25,.75)) 
bayes_R2(mod.review.Temp,probs = c(.25,.75))

save.image("..//forKNB/amer_plums_KNB.Rda")
