
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

setwd("~/Documents/git/proterant/FLOBUDS")
load(.Rdata)
d<-read.csv("first.event.dat.csv",header=TRUE)
###get rid of useless colummns
colnames(d)
d<-dplyr::select(d, -X.1)
d<-dplyr::select(d, -X)

#### make light differsent so variables  appear in order
d$Light<-ifelse(d$Light=="L","xL","S")
###name treatments numeric/continuous
d$photoperiod<-ifelse(d$Light=="xL",12,8)
d$temp_day<-ifelse(d$Force=="W",24,18)
d$temp_night<-ifelse(d$Force=="W",18,12)
d$chilldays<-ifelse(d$Chill==0,28,56)

###center predictors
d$p_cent<-d$photoperiod-mean(d$photoperiod)
d$f_cent<-d$temp_day-mean(d$temp_day)
d$c_cent<-d$chilldays-mean(d$chilldays)

###what is the distrbution of response variables
ggplot(d,aes(flo_day))+geom_density()
ggplot(d,aes(leaf_day))+geom_density()
ggplot(d,aes(Lbb_day))+geom_density()
ggplot(d,aes(Lexpand_day))+geom_density()
d<-unite(d, treatment, Force,Light, Chill, sep= "_",remove = FALSE)

###Basic plots for each phenophase#############################################
phased<-gather(d,phase,DOY,8:12)

bbview<-filter(phased, phase %in% c("Lbb_day","flo_day"))
expview<-filter(phased, phase %in% c("Lexpand_day","flo_day"))
leafview<-filter(phased, phase %in% c("leaf_day","flo_day"))


##filter to species that have goodish flowering and BB
bigsp<-filter(bbview, GEN.SPA %in% c("ACE.PEN", "COM.PER","COR.COR","ILE.MUC", "PRU.PEN","VAC.COR"))
ggplot(bigsp,aes(treatment,as.numeric(DOY)))+geom_point(aes(color=phase),size=0.5)+stat_summary(fun.data = "mean_cl_boot",aes(color=phase))+ggtitle("budburst vs. flower")+facet_wrap(~GEN.SPA)

###expansion and flowering
bigsp<-filter(expview, GEN.SPA %in% c("ACE.PEN", "COM.PER","COR.COR","ILE.MUC", "PRU.PEN","VAC.COR"))
ggplot(bigsp,aes(treatment,as.numeric(DOY)))+geom_point(aes(color=phase),size=0.5)+stat_summary(fun.data = "mean_cl_boot",aes(color=phase))+ggtitle("expansion vs. flower")+facet_wrap(~GEN.SPA)

###leafout and flowering
bigsp<-filter(leafview, GEN.SPA %in% c("ACE.PEN", "COM.PER","COR.COR","ILE.MUC", "PRU.PEN","VAC.COR"))
ggplot(bigsp,aes(treatment,as.numeric(DOY)))+geom_point(aes(color=phase),size=0.5)+stat_summary(fun.data = "mean_cl_boot",aes(color=phase))+ggtitle("leafout vs. flower")+facet_wrap(~GEN.SPA)

###################survival analysis###########Kaplan-Meier########################
viv<-filter(d,Dead.alive %in% c("A","?"))
table(viv$treatment)

table(viv$treatment)
table(d$treatment) 

vivo.full<-gather(viv,phase,DOY,9:12)

###do it for bud burst
vivo<-filter(vivo.full, phase %in% c("Lbb_day","flo_day"))
vivo$DOY<-ifelse(is.na(vivo$DOY),120,vivo$DOY)
vivo$surv<-ifelse(vivo$DOY==120,1,0)

vivo$floposs<-ifelse(vivo$Flo.poss.=="N",0,1)

###This is only twigs where we deemed flowering to even bee possible
vivo2<-filter(vivo,floposs==1)

###vivo2 is main datasheet for buds and flower
###chayim 2 is for leaf out and flowering
unique(vivo.full$phase)
hayim<-dplyr::filter(vivo.full, phase %in% c("leaf_day","flo_day"))
hayim$DOY<-ifelse(is.na(hayim$DOY),120,hayim$DOY)
hayim$surv<-ifelse(hayim$DOY==120,1,0)
hayim$floposs<-ifelse(hayim$Flo.poss.=="N",0,1)
###This is only twigs where we deemed flowering to even bee possible
hayim2<-filter(hayim,floposs==1)

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
ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

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








#### this model(above) is going to treat everything that didn't flower or leaf out as if the did at the same time :\


prior2<-get_prior(DOY | cens(surv) ~ phase+Light:phase+Chill:phase+Force:phase,
                 data = vivo2, family = weibull) 

###this model expludes things that died and couldn't have flowered
m1b<- brm(DOY | cens(surv) ~ phase+Light:phase+Chill:phase+Force:phase+(1+phase+Light:phase+Chill:phase+Force:phase|GEN.SPA),
          data = vivo2, family = weibull,inits = "0",
          iter= 3000,
          warmup = 2000,
          prior = prior2) 
summary(m1b)

Q<-as.data.frame(coef(m1b))
colnames(Q)

R<-dplyr::select(Q,"GEN.SPA.Estimate.phaseLbb_day", "GEN.SPA.Estimate.phaseflo_day.LightxL","GEN.SPA.Estimate.phaseLbb_day.LightxL","GEN.SPA.Estimate.phaseflo_day.Chill","GEN.SPA.Estimate.phaseflo_day.ForceW","GEN.SPA.Estimate.phaseLbb_day.Chill","GEN.SPA.Estimate.phaseLbb_day.ForceW")
R<-rownames_to_column(R,"GEN.SPA")
colnames(R)
colnames(R)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Force:Flo","Chill:Leaf","Force:Leaf")
R<-gather(R,"predictor","effect",2:8)
ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)
X<-rownames_to_column(X,"GEN.SPA")
colnames(X)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
X<-gather(X,"predictor","CIlow",2:8)

###high
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)
XX<-rownames_to_column(XX,"GEN.SPA")
colnames(XX)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
XX<-gather(XX,"predictor","CIhigh",2:8)

XXX<-full_join(X,XX)

Z<-full_join(XXX,R)

pd<- position_dodgev(height = 1)
Z$class<-NA
Z$class<-ifelse(Z$predictor=="Light:Flo",0,2)
Z$class<-ifelse(Z$predictor=="Light:Leaf",0,Z$class)

Z$class<-ifelse(Z$predictor=="Chill:Leaf",1,Z$class)
Z$class<-ifelse(Z$predictor=="Chill:Flo",1,Z$class)

Z$class<-ifelse(Z$predictor=="Phase",4,Z$class)

ggplot(Z,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model M1b:")+facet_wrap(~GEN.SPA)

Z2<-filter(Z, GEN.SPA %in% c("ACE.PEN","COM.PER","COR.COR","ILE.MUC","PRU.PEN","VAC.COR"))
ggplot(Z2,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model M1b")+facet_wrap(~GEN.SPA)

##########binary flower leaf or no model on suggestion of Lizzie


####Full model again leaf out vs. flowering


prior3<-get_prior(DOY | cens(surv) ~ phase+Light:phase+Chill:phase+Force:phase,
                  data = hayim2, family = weibull) 
###this model is best survivial
###do I need main effects for interpretation?
m2b<- brm(DOY | cens(surv) ~ phase+Light:phase+Chill:phase+Force:phase+(1+phase+Light:phase+Chill:phase+Force:phase|GEN.SPA),
          data = hayim2, family = weibull,inits = "0",
          iter= 3000,
          warmup = 2000,
          prior = prior3) 
summary(m2b)

Q<-as.data.frame(coef(m2b))
colnames(Q)

R<-dplyr::select(Q,"GEN.SPA.Estimate.phaseleaf_day", "GEN.SPA.Estimate.phaseflo_day.LightxL","GEN.SPA.Estimate.phaseleaf_day.LightxL","GEN.SPA.Estimate.phaseflo_day.Chill","GEN.SPA.Estimate.phaseflo_day.ForceW","GEN.SPA.Estimate.phaseleaf_day.Chill","GEN.SPA.Estimate.phaseleaf_day.ForceW")
R<-rownames_to_column(R,"GEN.SPA")
colnames(R)
colnames(R)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Force:Flo","Chill:Leaf","Force:Leaf")
R<-gather(R,"predictor","effect",2:8)
ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)
X<-rownames_to_column(X,"GEN.SPA")
colnames(X)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
X<-gather(X,"predictor","CIlow",2:8)

###high
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)
XX<-rownames_to_column(XX,"GEN.SPA")
colnames(XX)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
XX<-gather(XX,"predictor","CIhigh",2:8)

XXX<-full_join(X,XX)

Z<-full_join(XXX,R)

pd<- position_dodgev(height = 1)
Z$class<-NA
Z$class<-ifelse(Z$predictor=="Light:Flo",0,2)
Z$class<-ifelse(Z$predictor=="Light:Leaf",0,Z$class)

Z$class<-ifelse(Z$predictor=="Chill:Leaf",1,Z$class)
Z$class<-ifelse(Z$predictor=="Chill:Flo",1,Z$class)

Z$class<-ifelse(Z$predictor=="Phase",4,Z$class)

ggplot(Z,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model M2b:")+facet_wrap(~GEN.SPA)

Z2<-filter(Z, GEN.SPA %in% c("ACE.PEN","COM.PER","COR.COR","ILE.MUC","PRU.PEN","VAC.COR"))
ggplot(Z2,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model M2b")+facet_wrap(~GEN.SPA)
summary(m2b)
###Attempt to transform the real values
goober<-dplyr::select(Q, contains(".Estimate"))
goober<-rownames_to_column(goober,"GEN.SPA")
ncol(goober)
colnames(goober)<-c("GEN.SPA","Intercept","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
goober<-gather(goober,"predictor","effect",2:9)
goober$adj_effect<-exp(goober$effect)

summary(m2b) ###WHen exponentiated the #s are very high. Let's see if we just do a leaf based model
viv<-filter(d,Dead.alive %in% c("A","?"))
table(viv$treatment)

table(viv$treatment)
table(d$treatment) 

vivo.full<-gather(viv,phase,DOY,9:12)

###do it for bud burst only######################################
vivoo<-filter(vivo.full, phase %in% c("Lbb_day"))
vivoo$DOY<-ifelse(is.na(vivoo$DOY),120,vivoo$DOY)
vivoo$surv<-ifelse(vivoo$DOY==120,1,0)
table(vivoo$DOY) #182 out of 501 are censored

priorz<-get_prior(DOY | cens(surv) ~ Light+Chill+Force,
                  data = vivoo, family = weibull) 

###this model only expludes things that died not couldn't ahve flowered
###do I need main effects for interpretation?
mLonly<- brm(DOY | cens(surv) ~ Light+Chill+Force+(1+Light+Chill+Force|GEN.SPA),
          data = vivoo, family = weibull,inits = "0",
          iter= 3000,
          warmup = 2000,
          prior = priorz) 
summary(mLonly)
exp(5.39)
######################
hayim2$bin<-ifelse(hayim2$DOY==120,0,1)
prior4<-get_prior(bin ~ phase+Light:phase+Chill:phase+Force:phase,
                  data = hayim2, family=bernoulli(link = "logit") ) 
?brmsfamily()
###this model is best survivial
###do I need main effects for interpretation?
mbinom<- brm(bin ~ phase+Light:phase+Chill:phase+Force:phase+(1+phase+Light:phase+Chill:phase+Force:phase|GEN.SPA),
          data = hayim2, family =bernoulli(link = "logit"),
          iter= 3000,
          warmup = 2000,
          prior = prior4)
summary(mbinom)
coef(mbinom)
Q<-as.data.frame(coef(mbinom))
colnames(Q)

R<-dplyr::select(Q,"GEN.SPA.Estimate.phaseleaf_day", "GEN.SPA.Estimate.phaseflo_day.LightxL","GEN.SPA.Estimate.phaseleaf_day.LightxL","GEN.SPA.Estimate.phaseflo_day.Chill","GEN.SPA.Estimate.phaseflo_day.ForceW","GEN.SPA.Estimate.phaseleaf_day.Chill","GEN.SPA.Estimate.phaseleaf_day.ForceW")
R<-rownames_to_column(R,"GEN.SPA")
colnames(R)
colnames(R)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Force:Flo","Chill:Leaf","Force:Leaf")
R<-gather(R,"predictor","effect",2:8)
ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)
X<-rownames_to_column(X,"GEN.SPA")
colnames(X)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
X<-gather(X,"predictor","CIlow",2:8)

###high
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)
XX<-rownames_to_column(XX,"GEN.SPA")
colnames(XX)<-c("GEN.SPA","Phase","Light:Flo","Light:Leaf","Chill:Flo","Chill:Leaf","Force:Flo","Force:Leaf")
XX<-gather(XX,"predictor","CIhigh",2:8)

XXX<-full_join(X,XX)

Z<-full_join(XXX,R)

pd<- position_dodgev(height = 1)
Z$class<-NA
Z$class<-ifelse(Z$predictor=="Light:Flo",0,2)
Z$class<-ifelse(Z$predictor=="Light:Leaf",0,Z$class)

Z$class<-ifelse(Z$predictor=="Chill:Leaf",1,Z$class)
Z$class<-ifelse(Z$predictor=="Chill:Flo",1,Z$class)

Z$class<-ifelse(Z$predictor=="Phase",4,Z$class)

ggplot(Z,aes(effect,predictor))+geom_point(aes(color=as.character(class)))+geom_errorbarh(aes(color=as.character(class),xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+ggtitle("Model binom")+facet_wrap(~GEN.SPA)

###questions: Equivlency between chilling and other treatments?
