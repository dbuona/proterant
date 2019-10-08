###this model explores changes in dichogamy in a subset of species

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(tidyverse)
library(brms)
library(rstan)
library(arm)
library(rstanarm)
library(tibble)
library(ggstance)

setwd("~/Documents/git/proterant/FLOBUDS")

d<-read.csv("dicogamy.csv",header=TRUE)
plot.raw.dat<-gather(d,sex,DOY,5:6)

pd=position_dodge2(width=.8,preserve="total")
ggplot(plot.raw.dat,aes(Force,DOY))+geom_point(aes(shape=Light,color=sex),size=0.8)+stat_summary(aes(shape=Light,color=sex),position=pd)+facet_grid(Chill~GEN.SPA)+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Forcing")+ggthemes::theme_base()



d<-unite(d,treatment,Force,Light,Chill,sep="",remove=FALSE)
d$Light<-ifelse(d$Light=="L",1,0)
d$Force<-ifelse(d$Force=="W",1,0)

###name treatments numeric/continuous if you wany
d$photoperiod<-ifelse(d$Light==1,12,8)
d$temp_day<-ifelse(d$Force==1,24,18)
d$temp_night<-ifelse(d$Force==1,18,12)
d$chilldays<-ifelse(d$Chill==0,28,56)


d$floM<-ifelse(is.na(d$flo_dayM),0,1)
d$floF<-ifelse(is.na(d$flo_dayF),0,1)
cor<-filter(d,GEN.SPA=="COR.COR")
com<-filter(d,GEN.SPA=="COM.PER")

### likelihood of flowering

com
ladylikelinood<-brm(floF~Force+Chill+Light+Force:Chill+Force:Light+Chill:Light,data=cor, family=bernoulli(link="logit"))
summary(ladylikelinood)
brolikelinood<-brm(floM~Force+Chill+Light+Force:Chill+Force:Light+Chill:Light,data=cor, family=bernoulli(link="logit"))



summary(brolikelinood)
extract_coefs<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.10,.25,.75,0.90))),"Predictor")
}
brm(flo_dayM~Force+Chill+Light+Force:Chill+Force:Light+Chill:Light,data=com,iter=5000,warmup=4000)

lady<-extract_coefs(ladylikelinood)
bro<-extract_coefs(brolikelinood)

lady$phase<-"female"
bro$phase<-"male"

bothy<-rbind(lady,bro)
bothy<-filter(bothy,Predictor!="Intercept")
bothy$Predictor[bothy$Predictor=="Force"]<-"main:Force"
bothy$Predictor[bothy$Predictor=="Light"]<-"main:Light"
bothy$Predictor[bothy$Predictor=="Chill"]<-"main:Chill"
bothy$Predictor[bothy$Predictor=="Chill:Light"]<-"int:Chill:Light"
bothy$Predictor[bothy$Predictor=="Light:Force"]<-"int:Light:Force"
bothy$Predictor[bothy$Predictor=="Chill:Force"]<-"int:Chill:Force"

pd2=position_dodgev(height=0.3)
ggplot(bothy,aes(Estimate,Predictor))+geom_point(aes(color=phase,shape=phase),position=pd2, size=4)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")
dev.off()

lady.fect<-brm(flo_dayF~Force+Chill+Light+Force:Chill+Force:Light+Chill:Light,data=cor)
summary(lady.fect)
bro.fect<-brm(flo_dayM~Force+Chill+Light+Force:Chill+Force:Light+Chill:Light,data=cor)
summary(bro.fect)

lady<-extract_coefs(lady.fect)
bro<-extract_coefs(bro.fect)

lady$phase<-"female"
bro$phase<-"male"

bothy<-rbind(lady,bro)
bothy %>%
  arrange(Estimate) %>%
  mutate(Predictor = factor(Predictor, levels=c("Force:Chill",  "Force:Light","Chill:Light", "Chill", "Force", "Light","Intercept"))) %>%
  ggplot(aes(Estimate,Predictor))+geom_point(aes(color=phase,shape=phase),size=3,position=pd2)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),size=.5,linetype="dashed",stat="identity",height=0,position=pd2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),size=.5,height=0,position=pd2)+geom_vline(aes(xintercept=0),color="black")+theme_bw()+ggtitle("Cornus cornuta")

pd2=position_dodgev(height=0.3)
ggplot(bothy,aes(Estimate,Predictor))+geom_point(aes(color=phase,shape=phase),position=pd2, size=4)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+ggtitle("C.cornuta")


lady.fect.2<-brm(flo_dayF~Force+Chill+Light+Force:Chill+Force:Light+Chill:Light,data=com)
summary(lady.fect.2)
bro.fect.2<-lm(flo_dayM~Force+Chill+Light+GEN.SPA+Force:GEN.SPA+Light:GEN.SPA+Chill:GEN.SPA,data=d)
summary(bro.fect.2)
lm(flo_dayM~Force+Chill+Light+Force:GEN.SPA+Light:GEN.SPA+Chill:GEN.SPA,data=com)

summary(bro.fect.2)

ladylikelinood.cor<-brm(gynlike~p.z+f.z+c.z+p.z:f.z+p.z:c.z+f.z:c.z,data=cor, family=bernoulli(link="logit"))
summary(ladylikelinood.cor)
brolikelinood.cor<-brm(anthlike~p.z+f.z+c.z,data=cor, family=bernoulli(link="logit"))
summary(brolikelinood.cor)


lady<-brm(flo_dayF~p.z+f.z+c.z,data=com)
bro<-brm(flo_dayM~~p.z+f.z,data=com) ### doesn't flower without chilling
summary(lady)
summary(bro)
lady2<-brm(flo_dayF~p.z+f.z+c.z,data=cor)
bro2<-brm(flo_dayM~~p.z+f.z+c.z,data=cor)


summary(lady2)
summary(bro2)
gynsums<- d %>% group_by(GEN.SPA,treatment) %>%  summarise(mean.gyn=mean(flo_dayF, na.rm=TRUE)) 
anthsums<- d %>% group_by(GEN.SPA,treatment)%>%  summarise(mean.anth=mean(flo_dayM,na.rm=TRUE))
gynsums2<- d %>% group_by(GEN.SPA,treatment) %>%  summarise(sd.gyn=sd(flo_dayF, na.rm=TRUE)) 
anthsums2<- d %>% group_by(GEN.SPA,treatment)%>%  summarise(sd.anth=sd(flo_dayM,na.rm=TRUE))
dicho<-left_join(gynsums,gynsums2)
dich2<-left_join(anthsums,anthsums2)
dicho<-left_join(dicho,dich2)
dicho$dicogamy<-dicho$mean.anth-dicho$mean.gyn

ggplot(dicho, aes(treatment, dicogamy))+geom_boxplot(aes(color=GEN.SPA))
dicho$dichogamybin<-ifelse(dicho$dicogamy<=0,1,0)
ggplot(dicho, aes(treatment,dichogamybin))+geom_point()+facet_wrap(~GEN.SPA)

summary(lady)

         