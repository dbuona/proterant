####to make potential prunus manuscript figures

#Note  I think there is something missing in how i extract the posteriors. iving the equivelent of ranef instead of coef
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

setwd("~/Documents/git/proterant/investment/Input")
load("plummy.Rda")
##read in cleaned data
d.flo<-read.csv("input_clean/FLS_clean.csv")

###############################
#1) Characterize species FLSs
#2) Test historic drought (macro-ecological patterns)
#3) Test plasticity

#1 I think this is the best model because it "controls" for day of year of obs (earlier observations more likely be to hysternathough)
mod.ord.scale<-brm(bbch.v.scale~doy+(doy|specificEpithet),data=d.flo,family=cumulative("logit"), warmup = 2500,iter=4000)
#mod.ord.short<-brm(bbch.short~doy+(doy|specificEpithet),data=d.flo,family=cumulative("logit"), warmup = 2500,iter=4000)


save.image("plummy.Rda")

#new.data<-data.frame(d.flo%>% group_by(specificEpithet)%>% summarise(doy=median(doy)))
#new.data<-data.frame(specificEpithet=unique(d.flo$specificEpithet),doy=rep(median(d.flo$doy),13))
#new.data<-data.frame(specificEpithet=unique(d.flo$specificEpithet))

new.data<-data.frame(quant=rep(c( "0%" , "25%",  "50%",  "75%" ,"100%"),13),d.flo%>% group_by(specificEpithet)%>% summarise(doy=quantile(doy)))


#predy<-fitted(mod.ord.short,newdata = new.data,probs = c(.25,.75))
predy<-fitted(mod.ord.scale,newdata = new.data,probs = c(.25,.75))
predy<-cbind(new.data,predy)
#predy<-fitted(mod.ord.nodoy,newdata = new.data,probs = c(.25,.75))


predy2<-predy %>%tidyr::gather("phase","likelihood",4:27)
predy2$species2<-predy2$specificEpithet
#predy2<-predy %>%tidyr::gather("phase","likelihood",1:16) # 4 short
predy.est<-filter(predy2,str_detect(phase, "^Estimate"))
predy.error<-filter(predy2,str_detect(phase, "^Q25"))
predy.error2<-filter(predy2,str_detect(phase, "^Q75"))

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

colnames(errorlow)[3]<-"Q25"
colnames(errorhigh)[3]<-"Q75"


result$Q25<-errorlow$Q25
result$Q75<-errorhigh$Q75
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

colnames(errorlow)[5]<-"Q25"
colnames(errorhigh)[5]<-"Q75"


result$Q25<-errorlow$Q25
result$Q75<-errorhigh$Q75
result
result1$goo<-paste(result1$species,result1$phase)
result$goo<-paste(result$species,result$phase)

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



jpeg("..//Plots/ord_grandmedian.jpeg", width=11, height=4,unit="in",res=300)
ggplot(data=result,aes(bbch,likelihood))+geom_point()+
  geom_ribbon(aes(x=int,ymin=0,ymax=likelihood),fill="lightgray",alpha=0.6)+facet_wrap(~species2,nrow=2)+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  ggtitle("Median of Dataset")
dev.off()

jpeg("..//Plots/ord_spmedian.jpeg", width=11, height=4,unit="in",res=300)
ggplot(data=result,aes(bbch,likelihood))+geom_point()+
  geom_ribbon(aes(x=int,ymin=0,ymax=likelihood),fill="lightgray",alpha=0.6)+facet_wrap(~species2,nrow=2)+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  ggtitle("Median for each species")
dev.off()

jpeg("..//Plots/ord_ginteronly.jpeg", width=11, height=4,unit="in",res=300)
ggplot(data=result,aes(bbch,likelihood))+geom_point()+
  geom_ribbon(aes(x=int,ymin=0,ymax=likelihood),fill="lightgray",alpha=0.6)+facet_wrap(~species2,nrow=2)+
  geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+ggthemes::theme_clean()+
  theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))+
  ggtitle("Intercept only")
dev.off()

result<-filter(result,quant!="100%")
#result<-filter(result,quant!="0%")



season<-as_labeller(c('0%'="Start of season",'25%'="Early season",'50%'="Mid season",'75%'="Late season"))
jpeg("..//Plots/ord_quants.jpeg", width=11, height=8,unit="in",res=300)
ggplot(data=result,aes(bbch,likelihood))+geom_point()+
  geom_ribbon(aes(x=int,ymin=0,ymax=likelihood,),alpha=0.3)+
  facet_grid(quant~species2,labeller=labeller(quant=season))+
  geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0))+ggthemes::theme_clean(base_size = 11)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))
dev.off()

#always1
sos<-c(1,1,1,1,1,1,1,1,1,1,0,0,1)
es<-c(0,0,1,1,0,1,1,0,0,0,0,0,1)
ms<-c(0,0,1,0,0,1,1,0,0,0,0,0,1)
ls<-c(0,0,0,0,0,0,1,0,0,0,0,0,1)


hystscore<-data.frame(specificEpithet=sort(unique(d.flo$specificEpithet)),sos,es,ms,ls)
hystscore$score <- rowSums(hystscore[2:5])

### now pdsi

#### can do this one for 
d.pdsi<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi.csv")
d.petal<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/petal_clean.csv")
d.fruit<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruitsize_clean.csv")
d.phen<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/fruit_phen.csv")
d.cold<-read.csv("~/Documents/git/proterant/investment/Input/input_clean/pruno_clean_pdsi_wint.csv")


#mmacropatters
d.pdsi<-left_join(d.pdsi,hystscore)
d.petal<-left_join(d.petal,hystscore)
d.fruit<-left_join(d.fruit,hystscore)
d.fruit<-filter(d.fruit,fruit_type=="fleshy")
d.phen<-left_join(d.phen,hystscore)
d.cold<-left_join(d.cold,hystscore)

#z.score everything for future analyses
zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
d.pdsi$pdsi.z<-zscore(d.pdsi$pdsi)
d.petal$petal.z<-zscore(d.petal$pental_lengh_mm)
d.fruit$fruit.z<-zscore(d.fruit$fruit_diam_mm)
d.phen$phen.z<-zscore(d.phen$doy)

pdsiraw<-dplyr::select(d.pdsi,lat,lon,pdsi.z)
#### run zscore modele for measurement error model
pdsi.mod.z<-brm(pdsi.z~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
petalmod.z<- brm(petal.z~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
fruitmod.z<- brm(fruit.z~(1|id)+(1|specificEpithet),data=d.fruit,warmup=2500,iter=4000)
phenmod.z<- brm(phen.z~(1|specificEpithet),data=d.phen,warmup=2500,iter=4000)

pdsiout<-dplyr::select(as.data.frame(coef(pdsi.mod.z)),1:2)
petalout<-dplyr::select(as.data.frame(coef(petalmod.z)$specificEpithet),1:2)
fruitout<-dplyr::select(as.data.frame(coef(fruitmod.z)$specificEpithet),1:2)
phenout<-dplyr::select(as.data.frame(coef(phenmod.z)),1:2)

colnames(pdsiout)<-c("pdsi_mean","pdsi_se")
colnames(petalout)<-c("petal_mean","petal_se")
colnames(fruitout)<-c("fruit_mean","fruit_se")
colnames(phenout)<-c("phen_mean","phen_se")

phenout$specificEpithet<-rownames(phenout)
fruitout$specificEpithet<-rownames(fruitout)
petalout$specificEpithet<-rownames(petalout)
pdsiout$specificEpithet<-rownames(pdsiout)

newdat<-left_join(fruitout,pdsiout)
newdat<-left_join(newdat,petalout)
newdat<-left_join(newdat,phenout)

d.flo2<-left_join(d.flo,newdat)
d.flo2$doy.z<-zscore(d.flo2$doy)


d.flo3<-dplyr::left_join(d.flo2,pdsiraw)
d.flo3<-d.flo3[!duplicated(d.flo3), ]

###this is the prefered model maybe depending on error

#mod.ord.wpreds1<-brm(bbch.v.scale~doy.z+pdsi.z+me(petal_mean,petal_se)+me(fruit_mean,fruit_se)+(doy.z|specificEpithet),data=d.flo3,family=cumulative("logit"), warmup = 2500,iter=4000)
mod.ord.wpreds2<-brm(bbch.v.scale~doy.z+me(pdsi_mean,pdsi_se)+me(petal_mean,petal_se)+me(fruit_mean,fruit_se)+me(phen_mean,phen_se)+(doy.z|specificEpithet),data=d.flo2,family=cumulative("logit"), warmup = 3500,iter=5000)

###fit below over night
#mod.ord.wpreds.woah<-brm(bbch.v.scale~doy.z+me(pdsi_mean,pdsi_se)+me(petal_mean,petal_se)+me(fruit_mean,fruit_se)+
   #       me(pdsi_mean,pdsi_se):me(petal_mean,petal_se)+me(pdsi_mean,pdsi_se):me(fruit_mean,fruit_se)+me(fruit_mean,fruit_se):me(petal_mean,petal_se)+
 #           (doy.z|specificEpithet),data=d.flo2,family=cumulative("logit"), warmup = 3000,iter=4000)

mod.ord.wpreds3<-brm(bbch.v.scale~doy.z+pdsi.z+me(petal_mean,petal_se)+me(fruit_mean,fruit_se)*me(phen_mean,phen_se)+(doy.z|specificEpithet),data=d.flo3,family=cumulative("logit"), warmup = 3000,iter=5000)

fixef(mod.ord.wpreds3,probs=c(.25,.75))

get_variables(mod.ord.wpreds2)
output2<-mod.ord.wpreds2 %>%
  spread_draws(bsp_mepdsi_meanpdsi_se,bsp_mepetal_meanpetal_se,bsp_mefruit_meanfruit_se)
colnames(output2)
output2 <-output2 %>% tidyr::gather("var","estimate",4:6)
library(bayesplot)

aa<-ggplot(output2,aes(y = var, x = estimate)) +
  tidybayes::stat_eye(fill="darkorchid1")+ggthemes::theme_few()+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_y_discrete(limits = c("bsp_mefruit_meanfruit_se","bsp_mepetal_meanpetal_se","bsp_mepdsi_meanpdsi_se"),labels=c("fruit size","petal size","pdsi"))+
  ylab("Variable")+xlab("Effect estimate")

jpeg("..//Plots/cerasus_mus.jpeg",width=12, height=7, units = "in",res=300)
ggpubr::ggarrange(aa,d, labels = c("a)",""),heights=c(.6,.4),nrow=2,ncol=1)
dev.off()



####non zscored models for comparative analysis using grouping factors
pdsi.mod<-brm(pdsi~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
minpdsi.mod<-brm(pdsi.min~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
petalmod<- brm(pental_lengh_mm~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
fruitlmod<- brm(fruit_diam_mm~(1|id)+(1|specificEpithet),data=d.fruit,warmup=2500,iter=4000)
phenlmod<- brm(doy~(1|specificEpithet),data=d.phen,warmup=2500,iter=4000)
cold.mod<-brm(wintert~(1|specificEpithet),data=d.cold,warmup=2500,iter=4000)

##extrac and group posteriors I think these are doing raneffs
goober2<-pdsi.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(goober2)

fixef(pdsi.mod)
flooby<-left_join(goober2,hystscore)
flooby<-flooby%>% group_by(score)


  mean_qi(Intercept,.width=0.5)
flooby<-filter(flooby,!is.na(flooby))





###peta
gooberp<-petalmod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 

flooby2<-left_join(gooberp,hystscore)
flooby2<-flooby2%>% group_by(score)

###fruit
gooberf<-fruitlmod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 

flooby3<-left_join(gooberf,hystscore)
flooby3<-flooby3%>% group_by(score)

gooberph<-phenlmod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(gooberph)


flooby4<-left_join(gooberph,hystscore)
flooby4<-flooby4%>% group_by(score)


goobermin<-minpdsi.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(goobermin)


floobymin<-left_join(goobermin,hystscore)
floobymin<-floobymin%>% group_by(score)%>%
  mean_qi(Intercept,.width=0.5)
floobymin<-filter(floobymin,!is.na(floobymin))


goobercold<-cold.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(goobercold)


floobycold<-left_join(goobercold,hystscore)
floobycold<-floobycold%>% group_by(score)%>%
  mean_qi(Intercept,.width=0.5)
floobycold<-filter(floobycold,!is.na(floobycold))

a<-ggplot(flooby,aes(score,Intercept))+stat_eye()+
  ylab("Mean 120 year PDSI at \n collection sites")+
  scale_x_continuous(name ="Flowering-first grouping",breaks=c(0,1,2,3,4),
                   labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)

#ggplot(floobymin,aes(score,Intercept))+geom_point(size=3)+geom_errorbar(aes(ymin=.lower,ymax=.upper,width=0))+
 # ylab("Min PDSI at collect sites")+
  #scale_x_continuous(name ="Flowering-first grouping",
   #                  labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_base(base_size = 11)


b<-ggplot(flooby2,aes(score,Intercept))+stat_eye()+ylab("petal length")+xlab("FLS group")+scale_x_continuous(name ="Flowering-first grouping",breaks=c(0,1,2,3,4),
                                                                                                         labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)
c<-ggplot(flooby3,aes(score,Intercept))+stat_eye()+ylab("fruit diameter")+xlab("FLS group")+scale_x_continuous(name ="Flowering-first grouping",breaks=c(0,1,2,3,4),
                                                                                                            labels=c("Never","At start of season","Through early season","Through mid season","Through late season"))+ggthemes::theme_few(base_size = 11)

d<-ggpubr::ggarrange(a,b,c,nrow=1,labels=c("b)","c)","d)"))

ggplot(flooby4,aes(score,Intercept))+stat_eye()+ylab("fruit phenology")+xlab("FLS group")+scale_x_continuous(name ="Flowering-first grouping",breaks=c(0,1,2,3,4),
                                                                                                             labels=c("Never","At start \nof season","Through \nearly season","Through \nmid season","Through \nlate season"))+ggthemes::theme_few(base_size = 11)

#ggplot(floobycold,aes(score,Intercept))+geom_point(size=3)+geom_errorbar(aes(ymin=.lower,ymax=.upper,width=0))+
  ylab("min winter T")
#e<-ggpubr::ggarrange(b,c,d, ncol=3,nrow=1)
#f<-ggpubr::ggarrange(a,e,ncol=1,nrow=2)

###plastic
d.um<-d.flo
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
pdsi.dater<-tidyr::gather(ext,"year","pdsi",1:119)
class(pdsi.dater$year)

pdsi.dater$year<-as.integer(pdsi.dater$year)

head(pdsi.dater)
head(pdsi.dater)
joiner<-dplyr::select(d.flo,specificEpithet,lat,lon)

pdsi.counter<-left_join(pdsi.dater,joiner)
pdsi.counter$dry<-ifelse(pdsi.counter$pdsi<-2,1,0)  

  
pdsi.count<-pdsi.counter %>%group_by(specificEpithet)%>% summarize(obs = n(),
                                                                   dry = sum(pdsi < -2,na.rm=TRUE)
  
                                                                  prop =dry / obs)

#brm(dry~(1|specificEpithet),data=pdsi.counter,family="bernoulli")

d.um<-left_join(d.um,pdsi.dater)


### a better way than mean


moda<-brm(bbch.v.scale~pdsi+doy+(pdsi+doy|specificEpithet),data=d.um,family=cumulative("logit"))
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
yaya<-moda%>%
  spread_draws(b_pdsi,r_specificEpithet[specificEpithet,pdsi])
tat<-filter()
ggplot(yaya,aes())

goo$effect<-ifelse(goo$species=="Main Effect","main","species")
q<-ggplot(goo,aes(Estimate2,species))+geom_point(aes(size=effect))+
  geom_errorbarh(aes(xmin=Q252,xmax=Q752,height=0))+scale_size_manual(values=c(4,2))+
    geom_vline(xintercept = 0, color="red")+ggthemes::theme_clean(base_size = 11) + theme(axis.text.y = element_text(face=ifelse(goo$species=="Main Effect","bold","italic")))+
  scale_y_discrete(name ="species", 
                   limits=c("alleghaniensis", "americana"    ,  "angustifolia"  , "gracilis"   ,    "hortulana"    , 
                           "maritima"  ,     "mexicana"    ,   "munsoniana"   ,  "nigra"     ,     "rivularis"   ,  
                            "subcordata" ,    "texana"     ,    "umbellata", "Main Effect"))+theme(legend.position = "none")+xlab("Drought effect estimate")+
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
tree<-read.tree("~/Documents/git/proterant/investment/Input/plum.tre")
is.ultrametric(tree)
 ## make ultrametric

names.intree<-tree$tip.label # names the names
namelist<-unique(d.flo$specificEpithet)

to.prune<-which(!names.intree%in%namelist) #prun the tree
pruned.by<-drop.tip(tree,to.prune)
plotTree(pruned.by)# this is the tree


###what are the tip labels in pruned phylogeny?

mytree.names<-pruned.by$tip.label # did i get them all
intersect(namelist,mytree.names) #yes


### Q1 DO we bother with phylogeny? Only 7 species, and as you can see, weak signal
meanfls<-d.flo %>% dplyr::group_by(specificEpithet)%>% dplyr::summarise(meanFLS=mean(bbch.v.scale),sdFLS=sd(bbch.v.scale))##Take mean FLS values
d.phylo<-filter(d.flo,specificEpithet %in% mytree.names)

d.phylo<-left_join(d.phylo,meanfls)
##quick phylo sig
meanfls4phylo<-filter(meanfls,specificEpithet %in% mytree.names)
meanfls4phylo<-left_join(meanfls4phylo,hystscore)
phylosig(pruned.by,meanfls4phylo$meanFLS,method="lambda",nsim = 100, test=TRUE) ### 0.1r
phylosig(pruned.by,meanfls4phylo$score,method="lambda",nsim = 100, test=TRUE) ### 0.1r

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("ggtree")
library(ggtree)
meanfls4phylo$tip.labels<-meanfls4phylo$specificEpithet
full_join(pruned.by,meanfls4phylo)
p<-ggtree(pruned.by1,branch.length="none")

jpeg("..//Plots/phylosig1", width=11, height=6,unit="in",res=300)
p %<+% meanfls4phylo+geom_tiplab(hjust=-.5)+geom_tippoint(aes(color=meanFLS),size=5)+ xlim(0, 5)+geom_cladelabel(node=3, label="lambda=0.169", 0,offset=-3)+scale_color_viridis_b()
dev.off()
jpeg("..//Plots/phylosig2", width=11, height=6,unit="in",res=300)

p %<+% meanfls4phylo+geom_tiplab(hjust=-.5)+geom_tippoint(aes(color=as.factor(score)),size=5)+ xlim(0, 5)+geom_cladelabel(node=3, label="lambda=7.02e-05", 0,offset=-3)
dev.off()

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

