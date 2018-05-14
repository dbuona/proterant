
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
library(survminer)

setwd("~/Documents/git/proterant/FLOBUDS")

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
d$p_cent<-d$photoperiod/mean(d$photoperiod)
d$f_cent<-d$temp_day/mean(d$temp_day)
d$c_cent<-d$chilldays/mean(d$chilldays)

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
vivo<-filter(d,Dead.alive %in% c("A","?"))
table(vivo$treatment)

table(vivo$treatment)
table(d$treatment) 

vivo<-gather(vivo,phase,DOY,9:12)

###do it for bud burst
vivo<-filter(vivo, phase %in% c("Lbb_day","flo_day"))
vivo$DOY<-ifelse(is.na(vivo$DOY),120,vivo$DOY)
vivo$surv<-ifelse(vivo$DOY==120,0,1)


surv_object<-Surv(time=vivo$DOY, event = vivo$surv,type="right")
s1<-survfit(surv_object ~ treatment, data = vivo)
summary(s1)
ggsurvplot(s1, data =vivo, fun = "event")
m1<-survreg(Surv(time=vivo$DOY ,event=vivo$surv)~phase+Chill+Light+Force+Chill:phase+phase:Force+Light:phase,data=vivo, dist="gaussian")
summary(m1)

##survival model in brms Not working from work computer but brms is wonky here.


m1a<- brm(DOY | cens(surv) ~ phase+Light+Chill+Force,
          data = vivo, family = gaussian,inits = "0") 
??isFALSE()   
?brm()




####survival analysis workish, similar to

))#######################################################################
###first models: All predictors as catagorical treatment variables
modflo<-brm((flo_day) ~ Light+Chill+Force+Light:Chill+Light:Force+Force:Chill+(Light +Chill+Force+Light:Chill+Light:Force+Force:Chill|p|GEN.SPA),
               data = d, family = gaussian, 
               iter= 3000,
               warmup = 2000,
               cores = 4)

summary(modflo) # 2 divergent
coef(modflo)

modleaf<-brm(leaf_day ~ Light+Chill+Force+Light:Chill+Light:Force+Force:Chill+(Light +Chill+Force+Light:Chill+Light:Force+Force:Chill|GEN.SPA),
            data = d, family = gaussian, 
            iter= 3000,
            warmup = 2000,
            cores = 4)

summary(modleaf) #2 divergent
#########Extract mean estimate
Q<-as.data.frame(coef(modleaf))
colnames(Q)

R<-dplyr::select(Q, "GEN.SPA.Estimate.LightxL","GEN.SPA.Estimate.Chill","GEN.SPA.Estimate.ForceW",
                              "GEN.SPA.Estimate.LightxL.Chill","GEN.SPA.Estimate.LightxL.ForceW","GEN.SPA.Estimate.Chill.ForceW")
R<-rownames_to_column(R,"GEN.SPA")
colnames(R)
colnames(R)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")
R<-gather(R,"predictor","effect",2:7)
#ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

############## Extract credible invervals
X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)
colnames(X)<-c("GEN.SPA.Estimate.LightxL","GEN.SPA.Estimate.Chill","GEN.SPA.Estimate.ForceW",
               "GEN.SPA.Estimate.LightxL.Chill","GEN.SPA.Estimate.LightxL.ForceW","GEN.SPA.Estimate.Chill.ForceW")

colnames(X)<-c("Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")

X<-rownames_to_column(X,"GEN.SPA")
X<-gather(X,"predictor","CIlow",2:7)

######now high CI
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)
colnames(XX)<-c("GEN.SPA.Estimate.LightxL","GEN.SPA.Estimate.Chill","GEN.SPA.Estimate.ForceW",
               "GEN.SPA.Estimate.LightxL.Chill","GEN.SPA.Estimate.LightxL.ForceW","GEN.SPA.Estimate.Chill.ForceW")
colnames(XX)<-c("Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")

XX<-rownames_to_column(XX,"GEN.SPA")
XX<-gather(XX,"predictor","CIhigh",2:7)

XXX<-full_join(X,XX)

Z<-full_join(XXX,R)
#ggplot(Z, aes(effect,predictor))+geom_point()+geom_point(size=2.5)+geom_segment(aes(y=predictor,yend=predictor,x=CIlow,xend=CIhigh))+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)
### Instad of faceting try ggstance pack pd function made by harold in code cat sent

pd<- position_dodgev(height = 1)
#ggplot(Z, aes(effect,predictor))+geom_point(position=pd,aes(color=GEN.SPA))+geom_errorbarh(aes(color=GEN.SPA,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))


#####Try with flowers
Q<-as.data.frame(coef(modflo))
colnames(Q)

R<-dplyr::select(Q, "GEN.SPA.Estimate.LightxL","GEN.SPA.Estimate.Chill","GEN.SPA.Estimate.ForceW",
          "GEN.SPA.Estimate.LightxL.Chill","GEN.SPA.Estimate.LightxL.ForceW","GEN.SPA.Estimate.Chill.ForceW")

R<-rownames_to_column(R,"GEN.SPA")
colnames(R)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")
R<-gather(R,"predictor","Feffect",2:7)
#ggplot(R, aes(Feffect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

############## Extract credible invervals
X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)
colnames(X)<-c("GEN.SPA.Estimate.LightxL","GEN.SPA.Estimate.Chill","GEN.SPA.Estimate.ForceW",
               "GEN.SPA.Estimate.LightxL.Chill","GEN.SPA.Estimate.LightxL.ForceW","GEN.SPA.Estimate.Chill.ForceW")


X<-rownames_to_column(X,"GEN.SPA")
colnames(X)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")
X<-gather(X,"predictor","FCIlow",2:7)

######now high CI
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)
colnames(XX)<-c("GEN.SPA.Estimate.LightxL","GEN.SPA.Estimate.Chill","GEN.SPA.Estimate.ForceW",
                "GEN.SPA.Estimate.LightxL.Chill","GEN.SPA.Estimate.LightxL.ForceW","GEN.SPA.Estimate.Chill.ForceW")


XX<-rownames_to_column(XX,"GEN.SPA")
colnames(XX)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")

XX<-gather(XX,"predictor","FCIhigh",2:7)

XXX<-full_join(X,XX)

Z2<-full_join(XXX,R)
#ggplot(Z2, aes(Feffect,predictor))+geom_point()+geom_point(size=2.5)+geom_segment(aes(y=predictor,yend=predictor,x=FCIlow,xend=FCIhigh))+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)
### Instad of faceting try ggstance pack pd function made by harold in code cat sent
colnames(Z2)<-c("GEN.SPA","predictor","CIlow","CIhigh","effect")
Z2$class<-"flower"
#ggplot(Z2, aes(effect,predictor))+geom_point(position=pd,aes(color=GEN.SPA))+geom_errorbarh(aes(color=GEN.SPA,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))
colnames(Z)
Z$class<-"leaves"
colnames(Z)

Z3<-rbind(Z,Z2)
colnames(Z3)

ggplot(Z3, aes(effect,predictor))+geom_point(position=pd,aes(color=class))+geom_errorbarh(aes(color=class,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+facet_wrap(~GEN.SPA)

#just species with lots of flowering: 
Zsub<-filter(Z3, GEN.SPA %in% c("ACE.PEN","COM.PER","ILE.MUC","COR.COR","PRU.PEN","VAC.COR"))
ggplot(Zsub, aes(effect,predictor))+geom_point(position=pd,aes(color=class))+geom_errorbarh(aes(color=class,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+facet_wrap(~GEN.SPA)


#####more quantitative model#######

modAL<-brm(leaf_day ~ photoperiod +chilldays+temp_day+photoperiod:chilldays+photoperiod:temp_day+temp_day:chilldays+(photoperiod +chilldays+temp_day+photoperiod:chilldays+photoperiod:temp_day+temp_day:chilldays|GEN.SPA),
                   data = d, family = gaussian, 
                  iter= 3500,
                  warmup = 3000,
                  cores = 4)
summary(modAL) #11 divergent transitions
coef(modAL)

modAF<-brm(flo_day ~ photoperiod +chilldays+temp_day+photoperiod:chilldays+photoperiod:temp_day+temp_day:chilldays+(photoperiod +chilldays+temp_day+photoperiod:chilldays+photoperiod:temp_day+temp_day:chilldays|GEN.SPA),
           data = d, family = gaussian, 
           iter= 4500,
           warmup = 4000,
           cores = 4)
coef(modAF) #42 diverge
pp_check(modAF)
pp_check (modAL)
###plot them as above:



mod.centF<-brm(flo_day ~ p_cent +c_cent+f_cent+p_cent:c_cent+p_cent:f_cent+f_cent:c_cent+(1+ p_cent +c_cent+f_cent+p_cent:c_cent+p_cent:f_cent+f_cent:c_cent|p|GEN.SPA),
          data = d, family = gaussian, 
          iter= 4000,
          warmup = 3500,
          cores = 4)
summary(mod.centF) ##3 divergent
mod.centL<-brm(leaf_day ~ p_cent +c_cent+f_cent+p_cent:c_cent+p_cent:f_cent+f_cent:c_cent+(1+ p_cent +c_cent+f_cent+p_cent:c_cent+p_cent:f_cent+f_cent:c_cent|p|GEN.SPA),
               data = d, family = gaussian, 
               iter= 4000,
               warmup = 3500,
               cores = 4)
summary(mod.centL) #9 divergent
bayes_R2(mod.centF)

######plot it
Q<-as.data.frame(coef(mod.centL))
colnames(Q)

R<-dplyr::select(Q, "GEN.SPA.Estimate.p_cent","GEN.SPA.Estimate.c_cent",
          "GEN.SPA.Estimate.f_cent","GEN.SPA.Estimate.p_cent.c_cent",
          "GEN.SPA.Estimate.p_cent.f_cent","GEN.SPA.Estimate.c_cent.f_cent")

R<-rownames_to_column(R,"GEN.SPA")
colnames(R)
colnames(R)<-c("GEN.SPA","Main:Photo","Main:Force","Main:Chill","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")
R<-gather(R,"predictor","effect",2:7)
ggplot(R, aes(effect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

############## Extract credible invervals
X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)

colnames(X)
colnames(X)<-c("Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")

X<-rownames_to_column(X,"GEN.SPA")
X<-gather(X,"predictor","CIlow",2:7)

######now high CI
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)

colnames(XX)<-c("Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")

XX<-rownames_to_column(XX,"GEN.SPA")
XX<-gather(XX,"predictor","CIhigh",2:7)

XXX<-full_join(X,XX)

Z<-full_join(XXX,R)
#### Do it for flowers
Q<-as.data.frame(coef(mod.centF))
colnames(Q)

R<-dplyr::select(Q, "GEN.SPA.Estimate.p_cent","GEN.SPA.Estimate.c_cent",
                 "GEN.SPA.Estimate.f_cent","GEN.SPA.Estimate.p_cent.c_cent",
                 "GEN.SPA.Estimate.p_cent.f_cent","GEN.SPA.Estimate.c_cent.f_cent")

R<-rownames_to_column(R,"GEN.SPA")
colnames(R)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")
R<-gather(R,"predictor","Feffect",2:7)
#ggplot(R, aes(Feffect,predictor))+geom_point()+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)

############## Extract credible invervals
X<-dplyr::select(Q, contains(".2.5.ile"))
colnames(X)
X<-dplyr::select(X, -GEN.SPA.2.5.ile.Intercept)
ncol(X)
colnames(X)<-c("GEN.SPA.Estimate.p_cent","GEN.SPA.Estimate.c_cent",
          "GEN.SPA.Estimate.f_cent","GEN.SPA.Estimate.p_cent.c_cent",
               "GEN.SPA.Estimate.p_cent.f_cent","GEN.SPA.Estimate.c_cent.f_cent")


X<-rownames_to_column(X,"GEN.SPA")
colnames(X)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")
X<-gather(X,"predictor","FCIlow",2:7)

######now high CI
XX<-dplyr::select(Q, contains(".97.5.ile"))
colnames(XX)
XX<-dplyr::select(XX, -GEN.SPA.97.5.ile.Intercept)
ncol(XX)
colnames(XX)<-c("GEN.SPA.Estimate.p_cent","GEN.SPA.Estimate.c_cent",
          "GEN.SPA.Estimate.f_cent","GEN.SPA.Estimate.p_cent.c_cent",
                "GEN.SPA.Estimate.p_cent.f_cent","GEN.SPA.Estimate.c_cent.f_cent")


XX<-rownames_to_column(XX,"GEN.SPA")
colnames(XX)<-c("GEN.SPA","Main:Photo","Main:Chill","Main:Force","Int:PhotoxChill","Int:PhotoxForce","Int:ChillxForce")

XX<-gather(XX,"predictor","FCIhigh",2:7)

XXX<-full_join(X,XX)

Z2<-full_join(XXX,R)
#ggplot(Z2, aes(Feffect,predictor))+geom_point()+geom_point(size=2.5)+geom_segment(aes(y=predictor,yend=predictor,x=FCIlow,xend=FCIhigh))+geom_vline(aes(xintercept=0,color="red"))+facet_wrap(~GEN.SPA)
### Instad of faceting try ggstance pack pd function made by harold in code cat sent
colnames(Z2)<-c("GEN.SPA","predictor","CIlow","CIhigh","effect")
Z2$class<-"flower"
#ggplot(Z2, aes(effect,predictor))+geom_point(position=pd,aes(color=GEN.SPA))+geom_errorbarh(aes(color=GEN.SPA,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))
colnames(Z)
Z$class<-"leaves"
colnames(Z)

Z3<-rbind(Z,Z2)
colnames(Z3)
ggplot(Z3, aes(effect,predictor))+geom_point(position=pd,aes(color=class))+geom_errorbarh(aes(color=class,xmin=(CIlow), xmax=(CIhigh)), position=pd, size=.5, height =0, width=0)+geom_vline(aes(xintercept=0))+facet_wrap(~GEN.SPA)

