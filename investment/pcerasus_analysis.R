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

#require(mapdata); require(maptools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/Documents/git/proterant/investment/Input")
load("mods_whatReviewerswant.Rda")

#load("pcerasus.Rda")
##read in cleaned data
d.flo<-read.csv("input_clean/FLS_clean.csv") ##data
tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre") ##tree



range(d.flo$year)
ggplot(d.flo,aes(year))+geom_histogram()
##### give tree branch lengths
#is.ultrametric(tree)
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

#d.flo$logFLS<-log(d.flo$bbch.v.scale) ## make FLS linear

#d.flo$rtFLS<-sqrt(d.flo$bbch.v.scale) 
#ggplot(d.flo,aes(rtFLS))+geom_density()




d.flo$YEAR.hin<-ifelse(d.flo$year<=1980,1980,d.flo$year)
d.flo$YEAR.hin<-d.flo$YEAR.hin-1980




d.flo.check <- d.flo %>%
  group_by(species) %>%
  filter(!(abs(doy - median(doy)) > 3*sd(doy)))
ggplot(d.flo,aes(doy))+geom_histogram(bins=20)+facet_wrap(~species)

ggplot(d.flo.check,aes(doy))+geom_histogram(bins=20)+facet_wrap(~species)

d.flo.check <- d.flo.check  %>% group_by(species) %>%  mutate(range = min(doy))

jpeg("..//Plots/whatReviwerswant/seasonal_distrbn.jpeg",height=7,width=7,units='in',res=300)
ggplot(d.flo.check,aes(doy))+geom_histogram(bins=15)+facet_wrap(~species,scales="free_x")+
  ggthemes::theme_few()+theme(strip.text = element_text(face="italic"))+xlab("day of year collected")
dev.off()
diffs<-filter(d.flo,!id %in% d.flo.check$id)

mod.ord.4review<-brm(bbch.v.scale~YEAR.hin+doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
mod.ord.4review.nodoy<-brm(bbch.v.scale~YEAR.hin+(1|species)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

mod.ord.4review.nooutlier<-brm(bbch.v.scale~YEAR.hin+doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d.flo.check,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
mod.ord.4review.nodoy.nooutlier<-brm(bbch.v.scale~YEAR.hin+(1|species)+(1|gr(specificEpithet, cov = A)),data=d.flo.check,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))


#mod.log.gause.4review<-brm(logFLS~year+doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=gaussian(), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))
#mod.rt.gause.4review<-brm(rtFLS~year+doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=gaussian(), warmup = 3000,iter=4000,control=list(adapt_delta=0.99))

#pp_check(mod.log.gause.4review) ### not great
#pp_check(mod.rt.gause.4review) ### not great neighter
#pp_check(mod.ord.4review) ## nails it

#summary(mod.ord.4review)


fixef(mod.ord.4review,probs = c(.055,.25,.75,.955))
fixef(mod.ord.4review.nodoy,probs = c(.055,.25,.75,.955))

fixef(mod.ord.4review.nooutlier,probs = c(.055,.25,.75,.955))
fixef(mod.ord.4review.nodoy.nooutlier,probs = c(.055,.25,.75,.955))

coef(mod.ord.4review.nodoy)
simpdat<-data.frame(species=unique(d.flo.check$specificEpithet),specificEpithet=unique(d.flo.check$specificEpithet))
simpdat$YEAR.hin<-20
prezsim<-fitted(mod.ord.4review.nodoy.nooutlier,newdata=simpdat,probs = c(.055,.955))
prezsim<-cbind(simpdat,prezsim)
prezest<- prezsim %>%
  dplyr::select(species,contains("Estimate"))
colnames(prezest)[2:7]<-c(1,2,3,4,5,6)
prezest$index.nodoy<-prezest$`1`+prezest$`2`+prezest$`3`
prezest$index2.nodoy<-prezest$`1`+prezest$`2`
prezest$index3.nodoy<-prezest$`1`
prezest<-tidyr::gather(prezest,"stage","probability",2:7)

prezlow<- prezsim %>%
  dplyr::select(species,contains("Q5"))
colnames(prezlow)[2:7]<-c(1,2,3,4,5,6)
prezlow<-tidyr::gather(prezlow,"stage","low",2:7)

prezhigh<- prezsim %>%
  dplyr::select(species,contains("Q9"))
colnames(prezhigh)[2:7]<-c(1,2,3,4,5,6)
prezhigh<-tidyr::gather(prezhigh,"stage","high",2:7)

#prezhigh<-dplyr::select(prezhigh,-stage)
#prezlow<-dplyr::select(prezlow,-stage)

prezest<-left_join(prezest,prezlow)
prezest<-left_join(prezest,prezhigh)
prezest<-distinct(prezest)
prezest$stage<-as.numeric(prezest$stage)

prezest$class<-ifelse(prezest$stage <4,"hysteranthous","non-hysteranthous")
indy<-dplyr::select(prezest,species,index.nodoy,index2.nodoy,index3.nodoy)
indy<-distinct(indy)
indy$index.nodoy<-round(indy$index.nodoy,2)

prezest$bbch<-NA
prezest$bbch[which(prezest$stage==1 )]<- "BBCH 0"
prezest$bbch[which(prezest$stage==2 )]<- "BBCH 09"
prezest$bbch[which(prezest$stage==3)]<- "BBCH 11"
prezest$bbch[which(prezest$stage==4)]<- "BBCH 15"
prezest$bbch[which(prezest$stage==5 )]<- "BBCH 17"
prezest$bbch[which(prezest$stage==6 )]<- "BBCH 19"


jpeg("..//Plots/whatReviwerswant/sps_preds_nodoy4supp.jpeg",height=7,width=7,units='in',res=300)
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
dev.off()






review_pred2 <- mod.ord.4review.nodoy %>% 
  epred_draws(newdata =simpdat,ndraws=1000,seed = 1234)
review_pred2$class<-ifelse(review_pred2$.category %in% c(1,2,3),"hysteranthous","non-hysteranthous")



ggplot(review_pred2,aes(.category,.epred))+geom_line(size=0.1,aes(group=.draw,color=class),alpha=0.5)+facet_wrap(~species)





season<-d.flo.check%>% dplyr::group_by(specificEpithet,species) %>% summarise(start=min(doy),end=max(doy))
#lat<-d.flo%>% dplyr::group_by(specificEpithet,species) %>% summarise(lat=mean(lat,na.rm=TRUE))
sps<-season$species

new.dat<-data.frame()
for (z in c(1:length(sps))){
  dat <-filter(season,species==sps[z])
  goo<-data.frame(doy=seq(dat$start,dat$end,by=1),species=rep(sps[z],length(dat$start)))
  new.dat<-rbind(goo,new.dat)
}



new.dat$YEAR.hin<-20
new.dat$specificEpithet<-new.dat$species


review_pred <- mod.ord.4review.nooutlier %>% 
  epred_draws(newdata =new.dat,ndraws=1000,seed = 1234)
review_pred$class<-ifelse(review_pred$.category %in% c(1,2,3),"hysteranthous","non-hysteranthous")
review_pred$class2<-ifelse(review_pred$.category %in% c(1,2),"hysteranthous","non-hysteranthous")


###example<-
examps<-filter(review_pred,species %in% c("umbellata","americana","subcordata"))

fil1<-examps %>%group_by(species)%>% summarise(doy=max(doy))
fil1$season<-"end of season"
fil2<-examps %>%group_by(species)%>% summarise(doy=min(doy))   
fil2$season<-"beginning of season"
fil3<-examps %>%group_by(species)%>% summarise(doy=median(doy))            
fil3$season<-"mid"

fillz<-rbind(fil1,fil2)
fillz<-left_join(fillz,examps)

fillz$bbch<-NA
fillz$bbch[which(fillz$.category==1 )]<- "BBCH 0"
fillz$bbch[which(fillz$.category==2 )]<- "BBCH 09"
fillz$bbch[which(fillz$.category==3)]<- "BBCH 11"
fillz$bbch[which(fillz$.category==4)]<- "BBCH 15"
fillz$bbch[which(fillz$.category==5 )]<- "BBCH 17"
fillz$bbch[which(fillz$.category==6 )]<- "BBCH 19"


pollywantaplot<-ggplot(fillz,aes(bbch,.epred))+geom_line(aes(group=.draw,color=class),size=.01,alpha=0.4)+stat_summary(aes(),color="lightgray",size=.75,shape=1)+facet_grid(species~season)+
  ggthemes::theme_few(base_size =11)+
  scale_color_viridis_d(direction = -1,begin=.6,end=.1)+xlab("")+
  theme(axis.text.x = element_text(angle = 270,vjust =.4,size = 7))
  
pollywantbplot<-ggplot(fillz,aes(class,.epred/1000))+ geom_col(aes(color=class,fill=class))+
  facet_grid(factor(species,levels = c("umbellata", "americana","subcordata"))~season)+xlab("")+
ggthemes::theme_few(base_size =11)+theme(legend.title = element_blank())+
  scale_color_viridis_d(direction = -1,begin=.6,end=.1)+theme(strip.text.y = element_text(face="italic"))+
  scale_fill_viridis_d(direction = -1,begin=.6,end=.1)+ylab("likelihood")

pdf("..//Plots/whatReviwerswant/conceptfig.pdf")
ggpubr::ggarrange(pollywantbplot,pollywantaplot,common.legend = TRUE,ncol=1,heights = c(.8,1))
dev.off()
ggpubr::ggarrange(pollywantaplot,pollywantbplot,common.legend = TRUE,ncol=1)


check<-review_pred %>% group_by(species,class,.draw,doy)%>%summarise(pred=sum(.epred))
check2<-review_pred %>% group_by(species,class2,.draw,doy)%>%summarise(pred=sum(.epred))

sumz<-check %>%group_by(species,class)%>% summarise(index=round(mean(pred),2))
sumz<-filter(sumz,class=="hysteranthous")

sumz2<-check2 %>%group_by(species,class2)%>% summarise(index=round(mean(pred),2))
sumz2<-filter(sumz2,class2=="hysteranthous")

sumz<-arrange(sumz,-index)
sumz2<-arrange(sumz2,-index)
sumz$species
sumz2$species

cofs<-as.data.frame(fixef(mod.ord.4review.nooutlier,probs = c(.055,.25,.75,.945)))
cofs<-tibble::rownames_to_column(cofs,"predictor")
cofs<-filter(cofs,predictor %in% c("doy","YEAR.hin"))
p2a<-ggplot(cofs,aes(Estimate,predictor))+geom_point(size=4)+
  geom_errorbarh(aes(xmin=`Q5.5`,xmax=`Q94.5`),height=0,size=.5)+
  geom_errorbarh(aes(xmin=`Q25`,xmax=`Q75`),height=0,size=1)+geom_vline(xintercept=0)+
  scale_y_discrete(labels=c("day of season","year"))+xlab("Estimated effect")+ggthemes::theme_few()

p2b<-ggplot(check,aes(doy,pred))+geom_line(aes(color=class,group=.draw),size=0.3)+facet_wrap(~factor(species,levels=c("mexicana","umbellata","angustifolia","maritima","gracilis","americana","munsoniana","alleghaniensis","nigra","hortulana","texana","rivularis","subcordata")))+
  ggthemes::theme_few(base_size =11)+
  scale_color_viridis_d(direction = -1,begin=.6,end=.1)+geom_label(data=sumz,x=150,y=.75,aes(label=index))+
  xlab("day of season")+
  ylab("hysteranthy likelihood")+
  theme(strip.text = element_text(face = "italic"))+theme(legend.position = "bottom",legend.title = element_blank())

jpeg("..//Plots/whatReviwerswant/sps_preds.jpeg",height=7,width=7,units='in',res=300)
ggpubr::ggarrange(p2a,p2b, ncol=1,heights=c(.3,.7),labels=c("a)","b)"))
dev.off()







check.h<-filter(check,class=="hysteranthous")
p2c<-ggplot(check.h,aes(doy,pred))+geom_line(aes(color=class,group=.draw),size=0.01)+facet_wrap(~factor(species,levels=c("mexicana","umbellata","angustifolia","maritima","gracilis","americana","munsoniana","alleghaniensis","nigra","hortulana","texana","rivularis","subcordata")))+
  ggthemes::theme_few(base_size =11)+
  scale_color_viridis_d(direction = -1,begin=.6,end=.1)+geom_label(data=sumz,x=150,y=.75,aes(label=index))+
  xlab("day of season")+
  ylab("hysteranthy likelihood")+
  theme(strip.text = element_text(face = "italic"))+theme(legend.position = "bottom",legend.title = element_blank())

jpeg("..//Plots/whatReviwerswant/sps_preds_altview.jpeg",height=7,width=7,units='in',res=300)
ggpubr::ggarrange(p2a,p2c, ncol=1,heights=c(.3,.7),labels=c("a)","b)"))
dev.off()

pdsi.dat<-read.csv("input_clean/pruno_clean_pdsi.csv")
pdsi.dat$species<-pdsi.dat$specificEpithet
pdsi.sum<-pdsi.dat %>% group_by(species) %>% summarise(mean.pdsi=mean(pdsi,na.rm=TRUE), sd.pdsi = sd(pdsi,na.rm=TRUE),n.pdsi = n(),se.pdsi = sd.pdsi / sqrt(n.pdsi))

petal.dat<-read.csv("input_clean/petal_clean.csv")
petal.dat$species<-petal.dat$specificEpithet
petal.sum<-petal.dat %>% group_by(species) %>% summarise(mean.petal=mean(pental_lengh_mm,na.rm=TRUE), sd.petal = sd(pental_lengh_mm,na.rm=TRUE),n.petal = n(),se.petal = sd.petal / sqrt(n.petal))
dep.sum<-left_join(petal.sum,pdsi.sum)
sumz<-left_join(sumz,dep.sum)

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

sumz$pdsi.z<-zscore(sumz$mean.pdsi)
sumz$petal.z<-zscore(sumz$mean.petal)
sumz<-left_join(sumz,indy)

cor(sumz$pdsi.z,sumz$petal.z)
cor(sumz$mean.pdsi,sumz$mean.petal)

brms::get_prior(bf(index ~ pdsi.z*petal.z, phi ~1),data = sumz,family=Beta())
bprior <- c(prior_string("student_t(3,0,.25)", class = "b"),
            prior_string("student_t(3,0,.25)", class = "Intercept"))
         
mod.review.wants<- brms::brm(
  brms::bf(index ~ pdsi.z*petal.z,
     phi ~1),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 5000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr") 

mod.review.wants.doy<- brm(
  brms::bf(index.nodoy ~ pdsi.z*petal.z,
     phi ~ 1),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 5000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr") 

dex<-sumz$index
stancode(mod.review.wants)

xtable(fixef(mod.review.wants,probs = c(.055,.25,.75,.945)))
xtable(fixef(mod.review.wants.doy,probs = c(.055,.25,.75,.945)))

conditional_effects(mod.review.wants.doy,prob=.89)

launch_shinystan(mod.review.wants)
library(brms)
pp_check(mod.review.wants,ndraws = 1000) ## at least its the right shape
pp_check(mod.review.wants.doy,nsamples = 100)

mod.review.wants.abit<- brm(
  bf(index ~ pdsi.z+petal.z,
     phi ~1),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")




compromise<-left_join(pdsi.dat,sumz)

mod.review.ish<- brm(
  bf(index ~ pdsi*mean.petal,
     phi ~pdsi*mean.petal),
  data = compromise,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr")

summary(lm(index~pdsi*mean.petal,data=compromise))

pp_check(mod.review.ish)

fixef(mod.review.ish)

conditional_effects(mod.review.ish)


cofs.2<-as.data.frame(fixef(mod.review.wants,probs = c(.055,.25,.75,.945)))
cofs.2<-tibble::rownames_to_column(cofs.2,"predictor")
cofs.2<-filter(cofs.2,!predictor %in% c("phi_Intercept","Intercept"))

fp<-ggplot(cofs.2,aes(x=Estimate,y=factor(predictor,level=c("pdsi.z:petal.z","pdsi.z","petal.z"))))+geom_point(size=4)+
  geom_errorbarh(aes(xmin=`Q5.5`,xmax=`Q94.5`),height=0,size=0.5)+
  geom_errorbarh(aes(xmin=`Q25`,xmax=`Q75`),height=0,size=1)+geom_vline(xintercept=0)+
  scale_y_discrete(name="predictor",labels=c("PDSI x petal length","PDSI","petal length"))+
  ggthemes::theme_few()+xlab("standardized effect size estimate")


p1<-plot(conditional_effects(mod.review.wants,prob=.89,surface = TRUE,method = c("fitted"),plot=FALSE))


pi<-p1[[1]]+ggthemes::theme_few()+ylim(0,1)+ylab("hysteranthy \nindex")+xlab("PDSI")
pi
pii<-p1[[2]]+ggthemes::theme_few()+ylim(0,1)+ylab("hysteranthy \nindex")+xlab("petal length")
pii
piii<-p1[[3]]+ggthemes::theme_few()+ylab("petal length")+xlab("PDSI")+scale_fill_discrete(name="hysteranthy \nindex",type = "viridis")+theme(legend.position = "right")
piii
p4<-ggpubr::ggarrange(pi,pii,nrow=1,common.legend = TRUE,legend="bottom",labels=c("b)","c)"))

p5<-ggpubr::ggarrange(p4,piii,labels=c("","d)"),widths=c(1,.5))

jpeg("..//Plots/whatReviwerswant/hypoth_preds.jpeg",height=7,width=7,units='in',res=200)
ggpubr::ggarrange(fp,p4,ncol=1,labels=c("a)",""),heights=c(.7,1))
dev.off()

p1<-plot(conditional_effects(mod.review.wants, "pdsi.z", ordinal = TRUE,prob = .5,plot=FALSE))

p3<-plot(conditional_effects(FNAordz.phylo2, "inflor.z", ordinal = TRUE,prob=.5,plot=FALSE))
conditions <- make_conditions(FNAordz.phylo2, "inflor.z")
p4<-plot(conditional_effects(FNAordz.phylo2, "meanpdsi.z",conditions=conditions,ordinal = TRUE,prob=.5,plot=FALSE))
range(FNA.small$inflor.z)
p1<-p1[[1]]+ggthemes::theme_few()+scale_y_discrete(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"))+xlab(" mean PDSI")

#p2<-p2[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p3<-p3[[1]]+ggthemes::theme_few()+ylab("")+xlab("inflorescence size")+theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())

p4<-p4[[1]]+ggthemes::theme_few()+xlab("inflorescence size")+scale_y_discrete(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"))

+scale_color_manual(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"),values=c("hotpink","orange","lightgreen","darkgreen"))+
  scale_fill_manual(name="FLS",labels=c("flowers before leaves","flowers before/with leaves","flowers with leaves","flowers after leaves"),values=c("hotpink","orange","lightgreen","darkgreen"))+xlab("fruit size")

potty<-ggpubr::ggarrange(p1,p3,common.legend = TRUE,ncol=2,legend="bottom",widths = c(.8,.5))


####old way for suppliment
compromise<-left_join(pdsi.dat,sumz)
mod.pdsi.phylo<-brm(pdsi~index+(1|specificEpithet)+(1|gr(species, cov = A)),data=compromise,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs                    
compromise2<-left_join(petal.dat,sumz)

mod.petal.phylo<-brm(pental_lengh_mm~index+(1|specificEpithet)+(1|gr(species, cov = A)),data=compromise2,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs                    

save.image("mods_whatReviewerswant.Rda")

lines.nophylo<-mod.pdsi.phylo%>%
  spread_draws(b_Intercept,  b_index )

a<-ggplot()+
  geom_jitter(data=compromise,aes(index,pdsi),color="black",fill="black",alpha=0.6,size=0.1,width = 0.1,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_index),alpha=0.01,color="skyblue3")+
  geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_index)),color="navy",size=2)+
  #geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean PDSI at collection sites")+
  ggthemes::theme_few(base_size = 11)
 

linespetal<-mod.petal.phylo%>%
  spread_draws(b_Intercept,  b_index )
b<-ggplot()+
  geom_jitter(data=compromise2,aes(index,pental_lengh_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.1,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linespetal,aes(intercept=b_Intercept,slope=b_index),alpha=0.01,color="skyblue3")+
  geom_abline(data=linespetal,aes(intercept=mean(b_Intercept),slope=mean(b_index)),color="navy",size=2)+
  #geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean petal length (mm")+
  ggthemes::theme_few(base_size = 11)



jpeg("..//Plots/dataplots_SUPP.jpeg", width=12, height=4,unit="in",res=200)
ggpubr::ggarrange(a,b,labels=c("a)","b)"),nrow=1)
dev.off()




library(xtable)
xtable(d.flo.check%>% group_by(species) %>% count())
indexcomps<-dplyr::select(sumz,species,index,index.nodoy)
xtable(indexcomps)
samps_sum<-dplyr::select(sumz,species,n.petal,n.pdsi)

fls.samps<-d.flo.check%>% group_by(species) %>% count()
colnames(fls.samps)[2]<-"n.FLS"
sampsize<-left_join(fls.samps,samps_sum)
xtable(sampsize)

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
#mod.ord.scale.phlyo<-brm(bbch.v.scale~doy+(doy|species)+(doy|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.95,max_treedepth=20)) ##


if(FALSE){
mod.ord.scale.phlyo.a<-brm(bbch.v.scale~doy+(doy|species)+(1|gr(specificEpithet, cov = A)),data=d.flo,data2=list(A = A),family=cumulative("logit"), warmup = 3000,iter=4000,control=list(adapt_delta=0.95,max_treedepth=15)) ##


fixef(mod.ord.scale.phlyo.a)
ranef(mod.ord.scale.phlyo.a)

##predict the ordinal
new.data<-data.frame(quant=rep(c( "0%" , "25%",  "50%",  "75%" ,"100%"),13),d.flo%>% dplyr::group_by(specificEpithet,species)%>% dplyr::summarise(doy=quantile(doy)))


season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))



predy2<-fitted(mod.ord.scale.phlyo.a,newdata = new.data,probs = c(.025,.25,.75,.975))
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
#jpeg("..//Plots/ord_quants_phylo.jpeg", width=11, height=11,unit="in",res=200)
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

library(ggtree)
jpeg("..//Plots/phylosig2.jpeg", width=4, height=4,unit="in",res=300)
p<-ggtree(pruned.by,layout = "roundrect")
p %<+% FLSindexB+geom_tiplab(hjust=-.2,align=TRUE,fontface="italic")+geom_tippoint(aes(color=as.factor(hystscoreB)),size=5,shape=15)+ xlim(0, 3)+scale_color_viridis_d(option="turbo",name="Hysteranthy",  labels=c("Never","At start of season","Through early season","Through mid season","Through late season"))
dev.off()

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

}
####pdsi model
d<-read.csv("input_clean/pdsi_spi.csv")
d<-dplyr::select(d,-X)
d$species<-d$specificEpithet
d<-left_join(d,sumz)

pdsi.sp.modeled<-brm(pdsi~(1|species),data=d,family=gaussian(), warmup = 3000,iter=4000)
pdsi.sp.coefs<-as.data.frame(coef(pdsi.sp.modeled))
pdsi.sp.coefs<-tibble::rownames_to_column(pdsi.sp.coefs, var = "species")


gooberino<-left_join(dep.sum,sumz)

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

gooberino$pdsi.z<-zscore(gooberino$mean.pdsi)
gooberino$petal.z<-zscore(gooberino$mean.petal)
gooberino$pdsisez<-zscore(gooberino$se.pdsi)
gooberino$petalsez<-zscore(gooberino$se.petal)


summary(lm(ratio~mean.pdsi,data=gooberino))
summary(lm(ratio~mean.petal,data=gooberino))
summary(lm(ratio~pdsi.z+petal.z,data=gooberino))
summary(lm(ratio~pdsi.z*petal.z,data=gooberino))

bform2 <- bf(ratio ~  mi(pdsi.z),data=gooberino) + bf(pdsi.z| mi(pdsisez),data=gooberino)




library(car)
vif(lm(ratio~mean.pdsi+mean.petal,data=gooberino))

cor(pdsi.sp.coefs$ratio,pdsi.sp.coefs$species.Estimate.Intercept,use = "pairwise.complete.obs")

pdsi.mod.review<- brm(
  bf(ratio ~ pdsi,
     phi ~pdsi),
  data = d,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 5000, warmup = 4000,
  cores = 4, seed = 1234,backend = "cmdstanr") 



if(FALSE){
mod.pdsi.phylo<-brm(pdsi~ratio+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs


d<-left_join(d,FLSindex)


mod.pdsi.phylo<-brm(pdsi~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
mod.pdsi.phyloB<-brm(pdsi~hystscoreB+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
mod.pdsi.phyloC<-brm(pdsi~hystscoreC+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 

chico<-d %>% group_by(specificEpithet) %>% summarise(mean(pdsi,na.rm=TRUE))

exreme<-filter(d, hystscoreB %in% c(0,4))



ggplot(exreme,aes(as.factor(hystscoreB),pdsi))+geom_boxplot()
ggplot(exreme,aes(as.factor(hystscoreB),))+geom_boxplot()
geom_point(aes(color=specificEpithet))
range(d$pdsi,na.rm=TRUE)
#mod.min.pdsi.phylo<-brm(pdsi.min~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##runs
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

}
###other covariates
d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
d.petal<-dplyr::select(d.petal,-X)
d.petal<-left_join(d.petal,sumz)

colnames(d.petal)
petal.mod.review<- brm(
  bf(ratio ~ pental_lengh_mm,
     phi ~pental_lengh_mm),
  data = d.petal,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr") 


range(d.petal$pental_lengh_mm,na.rm=TRUE)
petal.fake<-data.frame(pental_lengh_mm=seq(0.5,14,by=0.5))

petal_pred <- petal.mod.review %>% 
  epred_draws(newdata =petal.fake,ndraws=1000)

plotB<-ggplot(petal_pred , aes(x=pental_lengh_mm, y = .epred,group=.draw)) +
  geom_line(size=.01)+ggthemes::theme_few()+
  ylab("")+
  xlab("petal length (mm)")+ylim(.2,1)

ggpubr::ggarrange(plotA,plotB)



sumz<-left_join(sumz,dep.sum)

cor(sumz$mean.petal,sumz$mean.pdsi) #.537 too colinear?

both.review<- brm(
  bf(ratio ~ mean.pdsi|mi(se.pdsi)+mean.petal|mi(se.petal),
     phi ~mean.pdsi|mi(se.pdsi)+mean.petal|mi(se.petal)),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr") 

marginal_effects(both.review)
means.review<- brm(
  bf(ratio ~ mean.petal,
     phi ~mean.petal),
  data = sumz,
  family = Beta(),
  control=list(adapt_delta=.99),
  chains = 4, iter = 4000, warmup = 3000,
  cores = 4, seed = 1234,backend = "cmdstanr") 


bothfake=data.frame(mean.petal=c(3,4,5,6,7),mean.pdsi=c(-1,-.5,0,.5,1))


both_pred <- both.review %>% 
  epred_draws(newdata =bothfake)

ggplot(both_pred , aes(x=mean.petal, y = .epred,group=.draw)) +
  geom_line(size=.01)
ggplot(both_pred , aes(x=mean.pdsi, y = .epred,group=.draw)) +
  geom_line(size=.01)


if(FALSE){
d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")

d.fruit<-filter(d.fruit,fruit_type=="fleshy")


d.fruit<-dplyr::select(d.fruit,-X)

d.petal<-left_join(d.petal,FLSindex)
d.fruit<-left_join(d.fruit,FLSindex)
d.fruit$species<-d.fruit$specificEpithet
}

d.petal$species<-d.petal$specificEpithet






mod.petal.phylo<-brm(pental_lengh_mm~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##tree depth
mod.petal.phyloB<-brm(pental_lengh_mm~hystscoreB+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) ##tree depth
mod.petal.phyloC<-brm(pental_lengh_mm~hystscoreC+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.petal,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99))

tabfl<-data.frame(t(round(fixef(mod.petal.phyloB,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfl2<-data.frame(t(round(fixef(mod.petal.phylo,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfl3<-data.frame(t(round(fixef(mod.petal.phyloC,prob=c(.025,.25,.75,.975))[2,],digits=3)))

tabfl<-rbind(tabfl,tabfl2,tabfl3)
tabfl$classification<-c("main analaysis","alternate 1","alternate 2")
tabfl$Hystanthous_if<-c("50% fl. likelihood  with BBCH 0 & 09","25% fl. likelihood with BBCH 0","40% fl. likelihood with BBCH 0 & 09")
tabfl$mod_variable<-"petal length"


fixef(mod.petal.phylo,probs = c(.25,.75))
fixef(mod.petal.phyloB,probs = c(.25,.75))
fixef(mod.petal.phyloC,probs = c(.25,.75))

mod.fruit.phylo<-brm(fruit_diam_mm~hystscoreA+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 
mod.fruit.phyloB<-brm(fruit_diam_mm~hystscoreB+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 
mod.fruit.phyloC<-brm(fruit_diam_mm~hystscoreC+(1|specificEpithet)+(1|gr(species, cov = A)),data=d.fruit,data2=list(A=A),family=gaussian(),warmup=3500,iter=4500,control=list(adapt_delta=0.99)) 

save.image("pcerasus.Rda")

tabfr<-data.frame(t(round(fixef(mod.fruit.phyloB,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfr2<-data.frame(t(round(fixef(mod.fruit.phylo,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfr3<-data.frame(t(round(fixef(mod.fruit.phyloC,prob=c(.025,.25,.75,.975))[2,],digits=3)))
tabfr<-rbind(tabfr,tabfr2,tabfr3)
tabfr$classification<-c("main analaysis","alternate 1","alternate 2")
tabfr$Hystanthous_if<-c("50% fl. likelihood  with BBCH 0 & 09","25% fl. likelihood with BBCH 0","40% fl. likelihood with BBCH 0 & 09")
tabfr$mod_variable<-"fruit diameter"



fixef(mod.fruit.phylo,probs = c(.25,.75))
fixef(mod.fruit.phyloB,probs = c(.25,.75))
fixef(mod.fruit.phyloC,probs = c(.25,.75))
##B


suptab<-rbind(tab,tabfl,tabfr)
colnames(suptab)
suptab<-suptab[, c(9, 7, 8, 1,2,3,4,5,6)]
xtable::xtable(suptab)

###plot all that We're choosing scenario B as the best measure of hysteranthy
lines.nophylo<-mod.pdsi.phyloB%>%
  spread_draws(b_Intercept,  b_hystscoreB )

a<-ggplot()+
  geom_jitter(data=d,aes(hystscoreB,pdsi),color="black",fill="black",alpha=0.6,size=0.1,width = 0.25,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_hystscoreB),alpha=0.01,color="skyblue3")+
  geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_hystscoreB)),color="navy",size=2)+
  #geom_abline(data=lines.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=lines.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Mean PDSI at collection sites")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly \nseason","Through \nmid \nseason","Through \nlate \nseason"))+ggthemes::theme_few(base_size = 11)


linespetal<-mod.petal.phyloB%>%
  spread_draws(b_Intercept,  b_hystscoreB )

linesfruit<-mod.fruit.phyloB%>%
  spread_draws(b_Intercept,  b_hystscoreB )
b<-ggplot()+
  geom_jitter(data=d.fruit,aes(hystscoreB,fruit_diam_mm),color="black",fill="black",alpha=0.6,size=0.3,width = 0.25,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linesfruit,aes(intercept=b_Intercept,slope=b_hystscoreB),alpha=0.01,color="skyblue3")+
  geom_abline(data=linesfruit,aes(intercept=mean(b_Intercept),slope=mean(b_hystscoreB)),color="navy",size=2)+
  #geom_abline(data=linesfruit.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #geom_abline(data=linesfruit.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Fruit diameter (mm)")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly \nseason","Through \nmid \nseason","Through \nlate \nseason"))+ggthemes::theme_few(base_size = 11)


c<-ggplot()+
  geom_jitter(data=d.petal,aes(hystscoreB,pental_lengh_mm),color="black",fill="black",alpha=0.6,size=0.1,width = 0.25,height=0)+
  #stat_eye(data=d,aes(score,pdsi),alpha=0.6,fill="grey50")+
  geom_abline(data=linespetal,aes(intercept=b_Intercept,slope=b_hystscoreB),alpha=0.01,color="skyblue3")+
  geom_abline(data=linespetal,aes(intercept=mean(b_Intercept),slope=mean(b_hystscoreB)),color="navy",size=2)+
  # geom_abline(data=linespetal.nophylo,aes(intercept=b_Intercept,slope=b_score),alpha=0.004,color="firebrick1")+
  #  geom_abline(data=linespetal.nophylo,aes(intercept=mean(b_Intercept),slope=mean(b_score)),color="firebrick1",size=2)+
  ylab("Petal Length (mm)")+
  scale_x_continuous(name ="Hysteranthy",breaks=c(0,1,2,3,4),
                     labels=c("Never","At start \nof season","Through \nearly \nseason","Through \nmid \nseason","Through \nlate \nseason"))+ggthemes::theme_few(base_size = 11)




jpeg("..//Plots/dataplots.jpeg", width=12, height=4,unit="in",res=200)
ggpubr::ggarrange(a,b,c,labels=c("a)","b)","c)"),nrow=1)
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


