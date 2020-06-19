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
d<-read.csv("budburst_survival_data.csv")
d$phasebin<-ifelse(d$phase=="Lbb_day.9.",0,1) ## flowering i 1

###z.score to compare with effect interaction of binary phase Model 3
d$p_zz<-d$photoperiod-mean(d$photoperiod)/(2*(sd(d$photoperiod)))
d$f_zz<-d$temp_day-mean(d$temp_day)/(2*(sd(d$temp_day)))
d$c_zz<-d$chilldays-mean(d$chilldays)/(2*(sd(d$chilldays)))

###make predictors binary
d$Light<-ifelse(d$Light=="S",0,1)
d$Force<-ifelse(d$Force=="C",0,1)

d<- transform(d,taxa_num=as.numeric(factor(GEN.SPA)))

###new dataset
d.nosurv<-filter(d,surv==0) ###remove no bursts
d.flo<-filter(d.nosurv,phase=="flo_day.60.")### jsut flowers
d.leaf<-filter(d.nosurv,phase=="Lbb_day.9.") #### jsut leaves

#######################models 0.A and 0.B##################################### many divergent trans not actively in use.
#data.list<-with(d.nosurv,
         #       list(y=DOY,
        #             sp=taxa_num,
       #             chill=Chill,
      #              force=Force,
     #               photo=Light,
  ###                n_sp=length(unique(d.flo$taxa_num))))


#mod.flo.1 = stan('stan/winter_2level_floorleaf.stan', data = data.list,  
 #           iter = 3000, warmup=2200)

#mod.leaf.noint = stan('nointer_2level_flowleaf.stan', data = data.list,  
 #                 iter = 6000, warmup=5000) 



####plot the raw data
d.figure<-filter(d.nosurv !GEN.SPA %in% c("ACE.RUB"))
ggplot(d.nosurv,aes(treatment,DOY,color=phase))+stat_summary()+facet_wrap(~GEN.SPA)+theme_base()



##########################Model 1 a and b, model flowering and buds seperate;y############################################

prerleaf2<-get_prior(DOY ~Chill+Light+Force,data = d.leaf, family = gaussian())
mod.leaf2<- brm(DOY ~ ~Chill+Light+Force+(1+~Chill+Light+Force|GEN.SPA),
               data = d.leaf, family = gaussian(),
               iter= 4000,
               warmup = 3400,
               prior=prerleaf2) 
summary(mod.leaf2)
pp_check(mod.leaf2)


prerleaf2.int<-get_prior(DOY~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = d.leaf, family = gaussian())
mod.leaf2.int<-brm(DOY ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(1+Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                  data = d.leaf, family = gaussian(),
                  iter= 4000,
                  warmup = 3400)   
summary(mod.leaf2.int)
pp_check(mod.leaf2.int)


prerflo2<-get_prior(DOY ~Chill+Light+Force,data = d.flo, family = gaussian())
mod.flo2<-brm(DOY~Chill+Light+Force+(1+~Chill+Light+Force|GEN.SPA),
                data = d.flo, family = gaussian(),
                iter= 4000,
                warmup = 3400,
                prior=prerflo2)
  
###This is good

prerflo2.int<-get_prior(DOY~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = d.flo, family = gaussian())
mod.flo2.int<-brm(DOY ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(1+Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
      data = d.flo, family = gaussian(),
      iter= 4000,
      warmup = 3000)   

summary(mod.flo2.int)
pp_check(mod.flo2.int)


###functions
extract_coefs<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.10,.25,.75,0.90))),"Predictor")
}
extract_ranef<-function(x){dplyr::select(rownames_to_column(as.data.frame(ranef(x, summary=FALSE,probs=c(0.10,.25,.75,0.90))),"GEN.SPA"),-c(3,9,15,21,27,33,39))
}
reduce_ranef<-function(x){dplyr::select(x, )
}

###############################################

flowy<-extract_coefs(mod.flo2.int)
leafy<-extract_coefs(mod.leaf2.int)
leafy$phase<-"foliate"
flowy$phase<-"floral"

bothy<-rbind(flowy,leafy)
bothy<-filter(bothy,Predictor!="Intercept")

bothy$Predictor[bothy$Predictor=="Chill:Light"]<-"xChill:Light"
bothy$Predictor[bothy$Predictor=="Light:Force"]<-"xLight:Force"
bothy$Predictor[bothy$Predictor=="Chill:Force"]<-"xChill:Force"

pd2=position_dodgev(height=0.4)
hm<-ggplot(bothy,aes(Estimate,Predictor))+geom_point(aes(color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgrey", "black"))

ran.flo<-extract_ranef(mod.flo2.int)
colnames(ran.flo)

a<-dplyr::select(ran.flo,1:6)
b<-dplyr::select(ran.flo,1,7:11)
c<-dplyr::select(ran.flo,1,12:16)
d<-dplyr::select(ran.flo,1,17:21)
e<-dplyr::select(ran.flo,1,22:26)
f<-dplyr::select(ran.flo,1,27:31)
g<-dplyr::select(ran.flo,1,32:36)

a$Predictor<-"xIntercept"
b$Predictor<-"chilling"
c$Predictor<-"photoperiod"
d$Predictor<-"forcing"
e$Predictor<-"int:chillxphoto"
f$Predictor<-"int:chillxforce"
g$Predictor<-"int:photoxforce"

call<-c("GEN.SPA","Estimate","Q10","Q25","Q75","Q90","Predictor")
colnames(a)<-call
colnames(b)<-call
colnames(c)<-call
colnames(d)<-call
colnames(e)<-call
colnames(f)<-call
colnames(g)<-call
flo.sps<-rbind(a,b,c,d,e,f,g)
flo.sps$phase<-"floral"
pd2=position_dodgev(height=0.6)
ggplot(flo.sps,aes(Estimate,Predictor))+geom_point(aes(color=GEN.SPA),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=GEN.SPA),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=GEN.SPA),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()

###now leaf estiamte
ran.leaf<-extract_ranef(mod.leaf2.int)

aa<-dplyr::select(ran.leaf,1:6)
bb<-dplyr::select(ran.leaf,1,7:11)
cc<-dplyr::select(ran.leaf,1,12:16)
dd<-dplyr::select(ran.leaf,1,17:21)
ee<-dplyr::select(ran.leaf,1,22:26)
ff<-dplyr::select(ran.leaf,1,27:31)
gg<-dplyr::select(ran.leaf,1,32:36)

aa$Predictor<-"xIntercept"
bb$Predictor<-"chilling"
cc$Predictor<-"photoperiod"
dd$Predictor<-"forcing"
ee$Predictor<-"int:chillxphoto"
ff$Predictor<-"int:chillxforce"
gg$Predictor<-"int:photoxforce"

call<-c("GEN.SPA","Estimate","Q10","Q25","Q75","Q90","Predictor")
colnames(aa)<-call
colnames(bb)<-call
colnames(cc)<-call
colnames(dd)<-call
colnames(ee)<-call
colnames(ff)<-call
colnames(gg)<-call
leaf.sps<-rbind(aa,bb,cc,dd,ee,ff,gg)
leaf.sps$phase<-"foliate"

sps.plot<-rbind(leaf.sps,flo.sps)
sps.plot<-gather(sps.plot,"ltype","lower",3:4)
sps.plot<-gather(sps.plot,"utype","upper",3:4)
sps.plot$interval<-NA
sps.plot$interval[sps.plot$ltype=="Q10"]<-"80"
sps.plot$interval[sps.plot$utype=="Q90"]<-"80"
sps.plot$interval[sps.plot$ltype=="Q25"]<-"50"
sps.plot$interval[sps.plot$utype=="Q75"]<-"50"
sps.plot<-filter(sps.plot,Predictor!="xIntercept")
pd2=position_dodgev(height=.9)
ggplot(sps.plot,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0),color="black")+theme_base() ##these need to be added to the global estimates



#######bernouli yes no flower leaf out #################################################################
d.flo.full<-filter(d,phase=="flo_day.60.")### jsut flowers
d.leaf.full<-filter(d,phase=="Lexpand_day.11.")

flo.yn<-get_prior(surv~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = d.flo.full, family = bernoulli(link="logit"))

mod.floyn.int<-brm(surv ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(1+Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                  data = d.flo.full, family = bernoulli(link="logit"),
                  iter= 4000,
                  warmup = 3000,
                  prior=flo.yn)   
summary(mod.floyn.int)

leaf.yn<-get_prior(surv~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = d.leaf.full, family = bernoulli(link="logit"))

mod.leafyn.int<-brm(surv ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(1+Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                   data = d.leaf.full, family = bernoulli(link="logit"),
                   iter= 4000,
                   warmup = 3000,
                   prior=leaf.yn)   
summary(mod.leafyn.int)

flowy2<-extract_coefs(mod.floyn.int)
leafy2<-extract_coefs(mod.leafyn.int)
leafy2$phase<-"foliate"
flowy2$phase<-"floral"

bothy2<-rbind(flowy2,leafy2)
bothy2<-filter(bothy2,Predictor!="Intercept")

bothy2$Predictor[bothy2$Predictor=="Chill:Light"]<-"xChill:Light"
bothy2$Predictor[bothy2$Predictor=="Light:Force"]<-"xLight:Force"
bothy2$Predictor[bothy2$Predictor=="Chill:Force"]<-"xChill:Force"
yn<-ggplot(bothy2,aes(Estimate,Predictor))+geom_point(aes(color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("hotpink", "black"))

gridExtra::grid.arrange(yn,hm,ncol=2)





leaf.yn.z<-get_prior(surv~p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z,data = d.leaf.full, family = bernoulli(link="logit"))

mod.leafyn.z<-brm(surv ~p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z+(1+p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z|GEN.SPA),
                    data = d.leaf.full, family = bernoulli(link="logit"),
                    iter= 4000,
                    warmup = 3000,
                    prior=leaf.yn.z)   
summary(mod.leafyn.z)
stop("Only here and above is relevant: everthing below is old code")
###################################model 3 ### combine flowering and leafing into a single model################################################################
table(d$GEN.SPA)


######survival model-- i don't really understand the output, and Im not include biologically to believe these
prior2<-get_prior(DOY | cens(surv) ~ phasebin+p_zz+c_zz+f_zz+p_zz:phasebin+c_zz:phasebin+f_zz:phasebin,
                  data = d, family = weibull()) 

mod.fe.z<- brm(DOY | cens(surv) ~ phasebin+p_zz+c_zz+f_zz+p_zz:phasebin+c_zz:phasebin+f_zz:phasebin+(1+phasebin+p_zz+c_zz+f_zz+p_zz:phasebin+c_zz:phasebin+f_zz:phasebin|GEN.SPA),
          data = d, family = weibull(),inits = "0",
          iter= 4500,
          warmup = 3800,
          prior = prior2) 
summary(mod.fe.z) ## 12 divergetn transitions

########not survival model
prior.gaus<-get_prior(DOY~phasebin+p_zz+c_z+f_zz+p_zz:phasebin+c_zz:phasebin+f_zz:phasebin, data = d.nosurv, family = gaussian())
mod.fe.z.nosurv<- brm(DOY ~ phasebin+p_zz+c_z+f_zz+p_zz:phasebin+c_zz:phasebin+f_zz:phasebin+(1+phasebin+p_zz+c_zz+f_zz+p_zz:phasebin+c_zz:phasebin+f_zz:phasebin|GEN.SPA),
               data = d.nosurv, family = gaussian(),
               iter= 4000,
               warmup =3400,
               prior=prior.gaus) 
summary(mod.fe.z.nosurv) ##6 divergent transitions
ranef(mod.fe.z.nosurv)

fixef(mod.fe.z.nosurv)
fixed<-rownames_to_column(as.data.frame(fixef(mod.fe.z.nosurv,probs=c(0.1,0.9,0.25,0.75))),"Parameter")
fixed<-filter(fixed,Parameter!="Intercept")
fixed<-filter(fixed,Parameter!="phasebin")
pd2=position_dodgev(height=0.2)
ggplot(fixed,aes(Estimate,Parameter))+geom_point(position=pd2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90),linetype="dashed",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgrey","black"))

d$bin<-ifelse(d$DOY=)





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



#for now, skip to line 103
prerleaf<-get_prior(DOY ~ p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z,data = d.leaf, family = gaussian())

#mod.leaf<- brm(DOY ~ p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z+(1+p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z|GEN.SPA),
#              data = d.leaf, family = gaussian(),
#             iter= 4000,
#            warmup = 3500,
#           prior=prerleaf)   ##83 divergetn

mod.leaf<- brm(DOY ~ p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z+(1|GEN.SPA),
               data = d.leaf, family = gaussian(),
               iter= 4000,
               warmup = 3500,
               prior=prerleaf)  

summary(mod.leaf)
pp_check(mod.leaf)

prerflo<-get_prior(DOY ~  p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z,data = d.flo, family = gaussian())
#mod.flo<- brm(DOY ~ p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z+(1+ p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z|GEN.SPA),
#              data = d.flo, family = gaussian(),
#             iter= 4000,
#             warmup = 3500,
#            prior=prerflo)    ## 41 diverget 

mod.flo<- brm(DOY ~ p_z+c_z+f_z+p_z:c_z+p_z:f_z+c_z:f_z+(1|GEN.SPA),
              data = d.flo, family = gaussian(),
              iter= 4000,
              warmup = 3500)    
summary(mod.flo)
pp_check(mod.flo)





