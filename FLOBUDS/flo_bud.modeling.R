##This is Dan's main analysis. 2018

####data leadin 6 Nov 2018
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(tidyverse)
library("brms")
library(rstan)
library(arm)
library(rstanarm)
library(tibble)
library(ggstance)
library(survival)
library(sur)
library(survminer)
library(ggthemes)
library("Hmisc")
setwd("~/Documents/git/proterant/FLOBUDS")
d<-read.csv("flo_exapand_survival_data.csv")

d.flo<-filter(d,phase=="flo_day.60.")
d.leaf<-filter(d,phase=="Lexpand_day.11.")

###leaf only
prior<-get_prior(DOY | cens(surv) ~ Light+Chill+Force,
                 data = d.leaf, family = weibull) 

mod.leaf<- brm(DOY | cens(surv) ~ Light+Chill+Force+(1+Light+Chill+Force|GEN.SPA),
               data = d.leaf, family = weibull,inits = "0",
               iter= 3000,
               warmup = 2000,
               prior = prior)    

summary(mod.leaf)








###model for flowering and leaf exansion
table(d$GEN.SPA)
prior2<-get_prior(DOY | cens(surv) ~ phase+p_z+c_z+f_z+p_z:phase+c_z:phase+f_z:phase,
                  data = d, family = lognormal()) 

mod.fe.z<- brm(DOY | cens(surv) ~ phase+p_z+c_z+f_z+p_z:phase+c_z:phase+f_z:phase+(1+phase+Light:phase+Chill:phase+Force:phase|GEN.SPA),
          data = d, family = lognormal(),inits = "0",
          iter= 3000,
          warmup = 2000,
          prior = prior2) 
summary(mod.fe.z)
comp.only<-dplyr::filter(d,GEN.SPA %in% c("COM.PER","COR.COR","ILE.MUC","VAC.COR","PRU.PEN","ACE.PEN"))

prior.comp<-get_prior(DOY | cens(surv) ~ phase+p_z+c_z+f_z+p_z:phase+c_z:phase+f_z:phase,
                                           data = comp.only, family = lognormal()) 
                         

mod.comp.z<- brm(DOY | cens(surv) ~ phase+p_z+c_z+f_z+p_z:phase+c_z:phase+f_z:phase+(1+phase+p_z:phase+c_z:phase+f_z:phase|GEN.SPA),
               data = comp.only, family = lognormal(),inits = "0",
               iter= 3000,
               warmup = 2000,
               prior = prior.comp)
summary(mod.comp.z)

comp<-as.data.frame(fixef(mod.comp.z))
comp<-rownames_to_column(comp,"predictor")
ggplot(comp,aes(Estimate,predictor))+geom_point(aes(color=GEN.SPA))+geom_errorbarh(aes(xmin=(Q2.5), xmax=(Q97.5)))+facet+wrap(~GEN.SPA)+geom_vline(aes(xintercept=0))


####try the model above for individual species




#plot this model
ggplot(Z,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model M1b:")+facet_wrap(~GEN.SPA)
###### This model has 13 divergent transistions and no real effects in the predictors.
#Try running it numerically
colnames(d)
prior<-get_prior(DOY | cens(surv) ~ phase+photoperiod+chilldays+temp_day+photoperiod:phase+chilldays:phase+temp_day:phase,
                  data = d, family = weibull) 

mod.fe2<- brm(DOY | cens(surv) ~ phase+photoperiod+chilldays+temp_day+photoperiod:phase+chilldays:phase+temp_day:phase+(1+phase+photoperiod+chilldays+temp_day+photoperiod:phase+chilldays:phase+temp_day:phase|GEN.SPA),
             data = d, family = weibull,inits = "0",
             iter= 3000,
             warmup = 2000,
             prior = prior) 
summary(mod.fe) ### 18 divergent transitiosn
coef(mod.fe)

#####flower and leaves seperately
prior<-get_prior(DOY | cens(surv) ~ photoperiod+chilldays+temp_day,
                 data = d.flo, family = weibull) 

mod.flo<- brm(DOY | cens(surv) ~ photoperiod+chilldays+temp_day+(1+photoperiod+chilldays+temp_day|GEN.SPA),
              data = d.flo, family = weibull,inits = "0",
              iter= 3000,
              warmup = 2000,
              prior = prior)
coef(mod.flo)

prior<-get_prior(DOY | cens(surv) ~ photoperiod+chilldays+temp_day,
                 data = d.leaf, family = weibull) 

mod.leaf<- brm(DOY | cens(surv) ~ photoperiod+chilldays+temp_day+(1+photoperiod+chilldays+temp_day|GEN.SPA),
              data = d.leaf, family = weibull,inits = "0",
              iter= 3000,
              warmup = 2000,
              prior = prior)              

L<-as.data.frame(coef(mod.leaf))
Fl<-as.data.frame(coef(mod.flo))

summary(mod.leaf)
summary(mod.flo)
L$class<-"leafout"
Fl$class<-"flowering"
output<-rbind(L,Fl)

###leaf output
leaf<-dplyr::select(L, contains("Estimate"))
leaf<-rownames_to_column(leaf,"GEN.SPA")
colnames(leaf)<-c("GEN.SPA","Intercept","photo","chill","force")
leaf<-gather(leaf,"predictor","estimate",2:5)

low<-dplyr::select(L, contains("Q2.5"))
low<-rownames_to_column(low,"GEN.SPA")
colnames(low)<-c("GEN.SPA","Intercept","photo","chill","force")
low<-gather(low,"predictor","low",2:5)

high<-dplyr::select(L, contains("Q97.5"))
high<-rownames_to_column(high,"GEN.SPA")
colnames(high)<-c("GEN.SPA","Intercept","photo","chill","force")
high<-gather(high,"predictor","high",2:5)

X<-left_join(leaf,low)
XX<-left_join(X,high)
XX$class<-"leaf"
###flower output
flower<-dplyr::select(Fl, contains("Estimate"))
flower<-rownames_to_column(flower,"GEN.SPA")
colnames(flower)<-c("GEN.SPA","Intercept","photo","chill","force")
flower<-gather(flower,"predictor","estimate",2:5)

low<-dplyr::select(Fl, contains("Q2.5"))
low<-rownames_to_column(low,"GEN.SPA")
colnames(low)<-c("GEN.SPA","Intercept","photo","chill","force")
low<-gather(low,"predictor","low",2:5)

high<-dplyr::select(Fl, contains("Q97.5"))
high<-rownames_to_column(high,"GEN.SPA")
colnames(high)<-c("GEN.SPA","Intercept","photo","chill","force")
high<-gather(high,"predictor","high",2:5)

Y<-left_join(flower,low)
YY<-left_join(Y,high)
YY$class<-"flower"

plotdat<-rbind(YY,XX)
plotdat<-filter(plotdat,predictor!="Intercept")
ggplot(plotdat,aes(estimate,predictor))+geom_point(aes(color=class))+geom_errorbarh(aes(xmin=(low), xmax=(high),color=class))+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0))


#######try a model with just COMPER
VACO<-filter(hayim2,GEN.SPA=="COM.PER")

VAC1<-brm(DOY | cens(surv) ~ phase+photoperiod+chilldays+temp_day+photoperiod:phase+chilldays:phase+temp_day:phase,
    data = VACO, family = weibull,inits = "0",
    iter= 3000,
    warmup = 2000) 
 summary(VAC1)
launch_shinystan(VAC1)




####ofset analysis
offset1<-spread(vivo2,phase,DOY)
offset1$offset<-offset1$Lbb_day-offset1$flo_day
offset1<-dplyr::filter(offset1,flo_day!=120)
offset1<-dplyr::filter(offset1,Lbb_day!=120)

prior.offset<-get_prior(offset~Light+Chill+Force,
                  data = offset1, family = gaussian()) 




###cant do a weibull of there is a zero response
m1.offset<- brm(offset ~Light+Chill+Force+(1+Light+Chill+Force|GEN.SPA),
          data = offset1, family = gaussian(),
          iter= 3000,
          warmup = 2000,
          prior = prior.offset) 
summary(m1.offset)
coef(m1.offset)
####plot it
Q<-as.data.frame(coef(m1.offset))
colnames(Q)

R<-dplyr::select(Q,c(1,5,9,13))
R<-rownames_to_column(R,"GEN.SPA")
colnames(R)
colnames(R)<-c("GEN.SPA","Intercept","Light","Chill","Force")
R<-gather(R,"predictor","effect",2:5)
ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+geom_errorbarh(aes(,xmin=(low), xmax=(high)))+facet_wrap(~GEN.SPA)

X<-dplyr::select(Q, contains("Q2.5"))
colnames(X)

X<-rownames_to_column(X,"GEN.SPA")
colnames(X)<-c("GEN.SPA","Intercept","Light","Chill","Force")
X<-gather(X,"predictor","CIlow",2:5)

###high
XX<-dplyr::select(Q, contains("Q97.5"))
colnames(XX)
XX<-rownames_to_column(XX,"GEN.SPA")
colnames(XX)<-c("GEN.SPA","Intercept","Light","Chill","Force")
XX<-gather(XX,"predictor","CIhigh",2:5)

XXX<-full_join(X,XX)

Z<-full_join(XXX,R)

pd<- position_dodgev(height = 1)
Z$class<-NA
Z$class<-ifelse(Z$predictor=="Light:Flo",0,2)
Z$class<-ifelse(Z$predictor=="Light:Leaf",0,Z$class)

Z$class<-ifelse(Z$predictor=="Chill:Leaf",1,Z$class)
Z$class<-ifelse(Z$predictor=="Chill:Flo",1,Z$class)

Z$class<-ifelse(Z$predictor=="Phase",4,Z$class)

ggplot(Z,aes(effect,predictor))+geom_point()+geom_errorbarh(aes(,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("OFFSET MODEL:")+facet_wrap(~GEN.SPA)

Z2<-filter(Z, GEN.SPA %in% c("ACE.PEN","COM.PER","COR.COR","ILE.MUC","PRU.PEN","VAC.COR"))
ggplot(Z2,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model M1b")+facet_wrap(~GEN.SPA)








