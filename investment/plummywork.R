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

setwd("~/Documents/git/proterant/investment/Input")

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
  geom_errorbar(aes(ymin=Q25,ymax=Q75,width=0))+ggthemes::theme_base(base_size = 11)+theme(axis.text.x = element_text(angle = 300,hjust=-0.1))+ theme(strip.text = element_text(face = "italic"))
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

#mmacropatters
d.pdsi<-left_join(d.pdsi,hystscore)
d.petal<-left_join(d.petal,hystscore)
d.fruit<-left_join(d.fruit,hystscore)
d.fruit<-filter(d.fruit,fruit_type=="fleshy")
d.phen<-left_join(d.phen,hystscore)

pdsi.mod<-brm(pdsi~(1|specificEpithet),data=d.pdsi,warmup=2500,iter=4000)
petalmod<- brm(pental_lengh_mm~(1|id)+(1|specificEpithet),data=d.petal,warmup=2500,iter=4000)
fruitlmod<- brm(fruit_diam_mm~(1|id)+(1|specificEpithet),data=d.fruit,warmup=2500,iter=4000)
phenlmod<- brm(doy~(1|specificEpithet),data=d.phen,warmup=2500,iter=4000)


##extrac and group posteriors I think these are doing raneffs
goober2<-pdsi.mod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(goober2)


flooby<-left_join(goober2,hystscore)
flooby<-flooby%>% group_by(score)%>%
  mean_qi(Intercept,.width=0.5)
flooby<-filter(flooby,!is.na(flooby))

###peta
gooberp<-petalmod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 

flooby2<-left_join(gooberp,hystscore)
flooby2<-flooby2%>% group_by(score)%>%
  mean_qi(Intercept,.width=0.5)
flooby2<-filter(flooby2,!is.na(flooby2))

###fruit
gooberf<-fruitlmod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 

flooby3<-left_join(gooberf,hystscore)
flooby3<-flooby3%>% group_by(score)%>%
  mean_qi(Intercept,.width=0.5)
flooby3<-filter(flooby3,!is.na(flooby3))

gooberph<-phenlmod%>%
  spread_draws(r_specificEpithet[specificEpithet,term])%>%
  tidyr::spread(term,r_specificEpithet) 
colnames(gooberph)


flooby4<-left_join(gooberph,hystscore)
flooby4<-flooby4%>% group_by(score)%>%
  mean_qi(Intercept,.width=0.5)
flooby4<-filter(flooby4,!is.na(flooby4))


a<-ggplot(flooby,aes(score,Intercept))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper,width=0))+ylab("mean pdsi")+xlab("FLS group")+ggthemes::theme_base(base_size = 11)
b<-ggplot(flooby2,aes(score,Intercept))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper,width=0))+ylab("petal length")+xlab("FLS group")+ggthemes::theme_base(base_size = 11)
c<-ggplot(flooby3,aes(score,Intercept))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper,width=0))+ylab("fruit diameter")+xlab("FLS group")+ggthemes::theme_base(base_size = 11)
d<-ggplot(flooby4,aes(score,Intercept))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper,width=0))+ylab("fruit phenology")+xlab("FLS group")+ggthemes::theme_base(base_size = 11)

ggpubr::ggarrange(a,b,c,d, ncol=2,nrow=2)


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
goo<-tidyr::gather(ext,"year","pdsi",1:119)
class(goo$year)

goo$year<-as.integer(goo$year)

d.um<-left_join(d.um,goo)

moda<-brm(bbch.v.scale~pdsi+doy+(pdsi+doy|specificEpithet),data=d.um,family=cumulative("logit"))
coef(moda,probs = c(.25,.75))
