
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

library(ggplot2)
library("brms")
library(rstan)

library(ggstance)
library(dplyr)
library(ggthemes)
library("Hmisc")
library(brms)
library(broom)
library(RColorBrewer)
library("tidybayes")

setwd("~/Documents/git/proterant/FLOBUDS")

dat<-read.csv("input/dicogamy.csv",header = TRUE)

dat$Light<-ifelse(dat$Light=="S",0,1)
dat$Force<-ifelse(dat$Force=="C",0,1)
dat$treatment<- paste(dat$Force,dat$Light,dat$Chill)

dat$staminate<-ifelse(is.na(dat$flo_dayM),0,1)
dat$pistilate<-ifelse(is.na(dat$flo_dayF),0,1)

dat$stamday<-dat$flo_dayM
dat$pistday<-dat$flo_dayF

dichogat<-dat
dichogat$dichogamy<-dichogat$stamday-dichogat$pistday
dichogat<-filter(dichogat,!is.na(dichogamy))
dichogatcomp<-filter(dichogat,GEN.SPA=="COM.PER")

table(dichogat$GEN.SPA)

dat<-tidyr::gather(dat,"phase","flower",14:15)
dat<-tidyr::gather(dat,"phase2","phenology",14:15)

#dichgoamy model

#dichogmod<-brm(phenology~phase*GEN.SPA+(1|name),data=dat)




#mod.likem<-brm(likeM~Chill+(1|name),data=comperL,family = "bernoulli")
#mod.likef<-brm(likeF~Chill+(1|name),data=comperL,family = "bernoulli")
comper<-filter(dat,GEN.SPA=="COM.PER")
corcor<-filter(dat,GEN.SPA!="COM.PER")


unique(comper$name)
###cor cor
##corphen<-brm(phenology~Chill*Force*Light*phase+(1|name),
  #                   warmup=4000,iter=5000,control=list(adapt_delta=0.95),
   #                  data=corcor,family = "gaussian") ##converges

#compphen<-brm(phenology~Chill*Force*Light*phase+(1|name),
    #         warmup=4000,iter=5000,control=list(adapt_delta=0.95),
     #        data=comper,family = "gaussian") ##converges

posy<-position_dodge(width=0.1)
ggplot(dat,aes(GEN.SPA,phenology))+stat_summary(aes(color=phase),position=posy)+facet_grid(Chill~Light~Force~GEN.SPA,scales="free")


###addsreon regularization priors
prior <- c(
  prior(normal(0, 10), class = "b"),      # stronger regularization
  prior(normal(0, 10), class = "Intercept")
)

#compVern<-brm(flower~Chill*
 #               phase+Force*phase+Light*phase+(1|name),
  #            warmup=4000,iter=5000,control=list(adapt_delta=0.99),
   #           data=comper,family = bernoulli(link = "logit"),prior = prior)

#conditional_effects(compVern,prob = .5,)


compVern<-brm(flower~Chill*phase*Force*Light+(1|name),
              warmup=4000,iter=5000,control=list(adapt_delta=0.99),
              data=comper,family = bernoulli(link = "logit"),prior = prior)

#compPhen<-brm(phenology~Chill*phase*Force*Light+(1|name),
#             warmup=4000,iter=5000,control=list(adapt_delta=0.99),
 #             data=comper,family = gaussian(),prior = prior)

#conditional_effects(compPhen,prob = .5)
nuds2<-data.frame(Chill=rep(c(0,1),2),Force=rep(c(0,1),each=2),Light=rep(c(0,1),each=4),phase=rep(c("pistilate","staminate"),each=8))

library(tidybayes)
compy<-epred_draws(compVern,newdata = nuds2,ndraws = 1000,re_formula = NA)
compyplot<- compy %>%
  group_by(phase,Chill) %>%
  mean_qi(.epred, .width = c(0.5,0.89))

pd<-position_dodge(width=0.8)
p1<-ggplot(compyplot,aes(as.factor(Chill),.epred))+
  geom_bar(stat="identity",aes(fill=phase),position="dodge",alpha=0.6,color="black")+
  geom_errorbar(aes(ymin=.lower,ymax=.upper,group=phase,linewidth =as.factor(.width)),width=0.0,position=pd)+
  scale_x_discrete(name="vernalization",labels = c("30 days","60 days"))+
  scale_fill_manual(name="",values=c("gray90","gray30"))+ggthemes::theme_few()+ylab("Likelihood of flowering")+scale_linewidth_manual(name="",values=c(1,.25))

round(compyplot$.epred,2)

library(tidybayes)
library(dplyr)
library(ggplot2)

library(tidybayes)
library(dplyr)
library(ggplot2)
library(tidyr)

p2dat<-compVern %>%
  tidy_draws() %>%
  select(starts_with("b_")) %>% 
  pivot_longer(everything(),
               names_to = "parameter",
               values_to = "value") 

p2dat<-dplyr::filter(p2dat,!parameter %in%
    c("b_Intercept","b_phasestaminate"))
unique(p2dat$parameter)

p2dat$phsse<-ifelse(p2dat$parameter %in% c("b_Light","b_Force","b_Chill","b_Chill:Force","b_Chill:Light","b_Force:Light","b_Chill:Force:Light"),"pistillate","staminate")

p2dat$cue<-ifelse(p2dat$parameter %in% c("b_Light","b_phasestaminate:Light"),"photoperiod","forcing")
p2dat$cue<-ifelse(p2dat$parameter %in% c("b_Chill","b_Chill:phasestaminate"),"vernalization",p2dat$cue)
p2dat$cue<-ifelse(p2dat$parameter %in% c("b_Chill:Force","b_Chill:phasestaminate:Force"),"vernalization:forcing",p2dat$cue)
p2dat$cue<-ifelse(p2dat$parameter %in% c("b_Chill:Light","b_Chill:phasestaminate:Light"),"vernalization:photoperiod",p2dat$cue)
p2dat$cue<-ifelse(p2dat$parameter %in% c("b_Force:Light","b_phasestaminate:Force:Light" ),"forcing:photoperiod",p2dat$cue)
p2dat$cue<-ifelse(p2dat$parameter %in% c("b_Chill:Force:Light", "b_Chill:phasestaminate:Force:Light"),"vernalization:forcing:photoperiod",p2dat$cue)

desired_order <- c(
  "vernalization:forcing:photoperiod",
  "forcing:photoperiod",
  "vernalization:photoperiod",
  "vernalization:forcing",
  "photoperiod",
  "forcing",
  "vernalization"
)

p2dat$cue <- factor(p2dat$cue, levels = desired_order)

pd2<-position_dodgev(height=.4)
p2<-  ggplot(p2dat, aes(x = value, y = cue)) +
  stat_pointinterval(aes(color = phsse,shape=phsse), .width = c(.5, .89),position = pd2) +
  ggthemes::theme_few() +
  labs(
    x = "Effect size (log-odds)",
    y = NULL
  ) +
  scale_color_manual(values=c("gray70","gray30")) +
  geom_vline(xintercept = 0, linetype = "dashed")+xlim(-25,25)

pdf("../figure/vernalize.pdf")
ggpubr::ggarrange(p1,p2,ncol=1,common.legend = TRUE,heights = c(1.5,1),labels = c("a)","b)"))
dev.off()

#compVern<-brm(flower~Chill*phase+(1|name),
#              warmup=4000,iter=5000,control=list(adapt_delta=0.95),
 #             data=comper,family = "bernoulli")




#corVern<-brm(flower~Chill*phase+(1|name),
 #             warmup=4000,iter=5000,control=list(adapt_delta=0.95),
  #           data=corcor,family = "bernoulli")

#corVernForce<-brm(flower~Chill*phase*Force+(1|name),
 #            warmup=4000,iter=5000,control=list(adapt_delta=0.95),
  #          data=corcor,family = "bernoulli")


nuds2<-data.frame(Chill=rep(c(0,1),2),Force=rep(c(0,1),each=2),Light=rep(c(0,1),each=4),phase=rep(c("pistilate","staminate"),each=8))

library(tidybayes)
compy<-epred_draws(compVern,newdata = nuds2,ndraws = 1000,re_formula = NA)
compyplot<- compy %>%
  group_by(phase, Force,Light,Chill) %>%
  mean_qi(.epred, .width = 0.5)




compyplot$vernalization<-ifelse(compyplot$Chill==0,"low","high")
compyplot$forcing<-ifelse(compyplot$Force==0,"low","high")
compyplot$photoperiod<-ifelse(compyplot$Light==0,"8 hours","12 hours")

pd=position_dodge(width=.9)
ggplot(compyplot,aes(vernalization,.epred))+
  
  geom_bar(stat="identity",aes(fill=phase,color=photoperiod,linetype=forcing),position="dodge",alpha=0.7)+

  geom_errorbar(aes(ymin=.lower,ymax = .upper,fill=phase,color=photoperiod,linetype=forcing),width=0.1,position=pd)+scale_color_manual(values=c("black","tan"))+
  scale_fill_viridis_d(begin = 0.3,end=0.9)+ggthemes::theme_few()+ylab("Likelihood of flowering")
  
  ylab("Calendar growing season") +
  coord_cartesian(xlim = c(-20, 20)) +
  ggthemes::theme_few() +
  xlab("")+scale_color_viridis_d()+scale_shape_manual(name="year",values=c(0,1,2))



comppred<-fitted(compVern,newdata = nuds2,re_formula = NA,probs = c(.25,.75))
comppred<-cbind(nuds2,comppred)



corpred<-fitted(corVern,newdata = nuds2,re_formula = NA,probs = c(.25,.75))
corpred<-cbind(nuds2,corpred)


ggplot(comppred,aes(as.factor(Chill),Estimate*100))+
  geom_bar(stat="identity",aes(fill=phase),position="dodge",alpha=0.7,color="black")+ylim(0,100)+
  geom_errorbar(aes(ymin=Q25*100,ymax=Q75*100,group=phase),width=.1,position=pd)+
  scale_x_discrete(name="vernalizarion level",labels = c("low","high"))+
  scale_fill_viridis_d(begin = 0.3,end=0.9)+ggthemes::theme_few()+ylab("Likelihood of flowering")

ggplot(corpred,aes(as.factor(Chill),Estimate*100))+
  geom_bar(stat="identity",aes(fill=phase),position="dodge",alpha=0.7,color="black")+ylim(0,100)+
  geom_errorbar(aes(ymin=Q25*100,ymax=Q75*100,group=phase),width=.1,position=pd)+
  scale_x_discrete(name="vernalization level",labels = c("low","high"))+
  scale_fill_viridis_d(begin = 0.3,end=0.9)+ggthemes::theme_few()+ylab("Likelihood of flowering")


ggpubr::ggarrange(plota,plotb,common.legend=TRUE,labels = c("a)", "b)"))





corphen<-brm(phenology~Chill*phase+(1|name),
                 warmup=4000,iter=5000,control=list(adapt_delta=0.95),
                data=corcor,family = "gaussian") ##converges

conditional_effects(corphen)

compphen<-filter(comper,!is.na(flo_dayF) &!is.na(flo_dayM))
compphen<-distinct(compphen)
compphen$dichogamy<-compphen$flo_dayM-compphen$flo_dayF


corphen<-filter(corcor,!is.na(flo_dayF) &!is.na(flo_dayM))
corphen<-distinct(corphen)
corphen$dichogamy<-corphen$flo_dayM-corphen$flo_dayF


compphen$cat_dich<-ifelse(compphen$dichogamy<=0,0,1)

total<-compphen %>% group_by(Force,Light) %>% count()

dichog<-compphen %>% group_by(Force,Light,cat_dich) %>% count()
dichog<-filter(dichog,cat_dich==1)
colnames(dichog)[4]<-"n_pg"
cog<-left_join(total,dichog)
cog$perc<-cog$n_pg/cog$n

compphen$Force<-as.factor(compphen$Force)
compphen$Light<-as.factor(compphen$Light)


table(compphen$cat_dich)
mod.phen.cat<-brm(cat_dich~Force*Light+(Force*Light|name),
              warmup=4000,iter=5000,control=list(adapt_delta=0.95),
              data=compphen,family="bernoulli")


conditional_effects(mod.phen.cat,prob=.5)

mod.phen<-brm(dichogamy~Force*Light+(Force*Light|name),
                     warmup=4000,iter=5000,control=list(adapt_delta=0.95),
                     data=compphen) ##converges
conditional_effects(mod.phen,prob=.5,re_formula = ~(Force*Light|name))


mod.phen2<-brm(phenology~Force*Light*phase+(Force*Light*phase|name),
              warmup=4000,iter=5000,control=list(adapt_delta=0.95),
              data=compphen)


nudsP<-data.frame(name=rep(unique(compphen$name),4),Light=rep(c(0,1),each=14),Force=rep(c(0,1),28))
nudsP<-distinct(nudsP)

phenopredo<-epred_draws(mod.phen,newdata = nudsP,ndraws=1000)
phenopredo$grouper<-paste(phenopredo$.draw,phenopredo$name)

plotb<-ggplot(phenopredo,aes(Force,.epred))+geom_line(aes(group=grouper,color=name),size=0.01,alpha=.5)+
  geom_smooth(aes(color=name,fill=name),method="lm",se = TRUE,size=.5)+
  #geom_smooth(method="lm",se = FALSE,color="skyblue",size=2)+
  stat_pointinterval(aes(color=name),.width = c(.5))+
  facet_wrap(~Light,ncol=2)+xlim(0,1)+
  coord_cartesian(ylim=c(-40,60))+ggthemes::theme_few()+
  geom_hline(yintercept=0,color="black")+scale_color_viridis_d()+scale_fill_viridis_d()+theme(legend.position = "none")


ggpubr::ggarrange(plota,plotb,ncol=1)
comperL$Force<-as.factor(comperL$Force)
comperL$Chill<-as.factor(comperL$Chill)
comperL$Light<-as.factor(comperL$Light)

comperS<-filter(comperL,phase=="staminate")
comperP<-filter(comperL,phase!="staminate")

###try modeling male and female times seperately
mod.like.pistil<-brm(flower~Chill*Force*Light+(1|name),
                     warmup=4000,iter=5000,control=list(adapt_delta=0.95),
                     data=comperP,family = "bernoulli") ##converges

mod.like.stamen<-brm(flower~Chill*Force*Light+(1|name),
                     warmup=4000,iter=5000,control=list(adapt_delta=0.95),
                     data=comperS,family = "bernoulli") ###can't fit

fixef(mod.like.pistil,probs = c(.25,.75,.055,.945))



mod.like.stamen2<-brm(flower~Force*Light+(1|name),
                     warmup=4000,iter=5000,control=list(adapt_delta=0.95),
                     data=comperSC,family = "bernoulli") ###can't fit

mod.like.chillonly<-brm(flower~Chill*phase+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=comperL,family = "bernoulli")
conditional_effects(mod.like.chillonly,prob=.75)



corcorL$Force<-as.factor(corcorL$Force)
corcorL$Chill<-as.factor(corcorL$Chill)
corcorL$Light<-as.factor(corcorL$Light)
mod.like.chillonly.corcor<-brm(flower~Chill*phase+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=corcorL,family = "bernoulli")

conditional_effects(mod.like.chillonly.corcor,prob=.5)





conditional_effects(mod.like.pistil)

comperSC<-filter(comperS,Chill==1)



mod.like.stamen2<-brm(flower~Force*Light+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=comperSC,family = "bernoulli")


nuds<-data.frame(Chill=rep(c(0,1),2),Force=rep(c(0,1),each=2),Light=rep(c(0,1),each=4))
nuds2<-data.frame(Chill=c(0,1))
nuds3<-data.frame(Light=rep(c(0,1),2),Force=rep(c(0,1),each=2))

p1<-fitted(mod.like.pistil,newdata = nuds,probs = c(.05,.95,.25,.75),re_formula = NA)
p2<-fitted(mod.like.stamen2,newdata = nuds3,probs = c(.05,.95,.25,.75),re_formula = NA)

p1<-cbind(nuds,p1)

p2<-cbind(nuds3,p2)

ggplot(p2,aes(Force,Estimate))+geom_point(size=4,aes(color=as.factor(Light)),position=pd)+geom_errorbar(aes(ymin=Q25,ymax=Q75,color=as.factor(Light)),width=0,position=pd)
  
geom_errorbar(aes(ymin=Q5,ymax=Q95,color=phase),width=0,linetype="dotted",position=pd)+
  

p2$Chill<-1

p1$phase<-"pistillate"
p2$phase<-"staminate"
pp<-rbind(p1,p2)
pp$chill<-ifelse(pp$Chill==0,"low chilling","high chilling")
pp$force<-ifelse(pp$Force==0,"low forcing","high forcing")
pp$photoperiod<-ifelse(pp$Light==0,"short photoperiod","long photoperiod")

pd<-position_dodge(width=0.2)
ggplot(pp,aes(chill,Estimate))+geom_point(size=4,aes(color=phase,shape=phase),position=pd)+geom_errorbar(aes(ymin=Q25,ymax=Q75,color=phase),width=0,position=pd)+
   geom_errorbar(aes(ymin=Q5,ymax=Q95,color=phase),width=0,linetype="dotted",position=pd)+
  facet_grid(force~photoperiod)+ggthemes::theme_few()+ylim(0,1)+scale_color_viridis_d()+ylab("flowering likelihood")

ggplot(p1,aes(as.factor(Chill),Estimate))+geom_point(size=3)+geom_errorbar(aes(ymin=Q25,ymax=Q75),width=0)+
  geom_errorbar(aes(ymin=Q5,ymax=Q95),width=0,linetype="dotted")+
  facet_grid(Force~Light)+ggthemes::theme_base()+ylim(0,1)

conditional_effects(mod.like.pistil,prob= .5)
conditional_effects(mod.like.stamen,prob= .5)
conditional_effects(mod.like.stamen2,prob= .5)


conditional_effects(mod.likeC,prob= .5)
mod.likeC<-brm(flower~Chill*phase+Force*phase+Light*phase+(1|name),warmup=4000,iter=5000,data=comperL,family = "bernoulli")


mod.likeC<-brm(flower~Chill*phase+(1|name),warmup=4000,iter=5000,data=comperL,family = "bernoulli")

mod.likeFP<-brm(flower~Force*phase*Light+(1|name),warmup=4000,iter=5000,data=comperL,family = "bernoulli")
mod.likeP<-brm(flower~Light*phase+(1|name),warmup=4000,iter=5000,data=comperL,family = "bernoulli")
labs<-c(30,60)
labs2<-c(18,24)
labs3<-c(8,12)

p11<-plot(conditional_effects(mod.likeC, "Chill:phase", categorical = FALSE,prob = .5,plot=FALSE))
p11<-p11[[1]]+ggthemes::theme_few()+scale_color_viridis_d()+scale_fill_viridis_d()+scale_x_continuous(breaks=c(0,1),labels = labs)+ylab("flowering probability")+xlab("Days of chilling")

p22<-plot(conditional_effects(mod.likeF, "Force:phase", categorical = FALSE,prob = .5,plot=FALSE))
p22<-p22[[1]]+ggthemes::theme_few()+scale_color_viridis_d()+scale_fill_viridis_d()+scale_x_continuous(breaks=c(0,1),labels = labs2)+ylab("")+xlab("Daytime temperature")+ylim(0,1)


p33<-plot(conditional_effects(mod.likeP, "Light:phase", categorical = FALSE,prob = .5,plot=FALSE))
p33<-p33[[1]]+ggthemes::theme_few()+scale_color_viridis_d()+scale_fill_viridis_d()+scale_x_continuous(breaks=c(0,1),labels = labs3)+ylab("")+xlab("Photoperiod")+ylim(0,1)

ggpubr::ggarrange(p11,p22,p33,ncol=3,common.legend = TRUE)


corcorL$Force<-as.factor(corcorL$Force)
corcorL$Chill<-as.factor(corcorL$Chill)
corcorL$Light<-as.factor(corcorL$Light)

mod.likeCc<-brm(flower~Light*phase*Chill*Force+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=corcorL,family = "bernoulli")
conditional_effects(mod.likeCc,prob=.5)


mod.likeCc<-brm(flower~Chill*phase+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=corcorL,family = "bernoulli")
mod.likeFc<-brm(flower~Force*phase+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=corcorL,family = "bernoulli")
mod.likePc<-brm(flower~Light*phase+(1|name),warmup=4000,iter=5000,control=list(adapt_delta=0.95),data=corcorL,family = "bernoulli")


pCc<-plot(conditional_effects(mod.likeCc, "Chill:phase", categorical = FALSE,prob = .5,plot=FALSE))
pCc<-pCc[[1]]+ggthemes::theme_few()+scale_color_viridis_d()+scale_fill_viridis_d()+scale_x_continuous(breaks=c(0,1),labels = labs)+ylab("flowering probability")+xlab("Days of chilling")+ylim(0,1)

pFc<-plot(conditional_effects(mod.likeFc, "Force:phase", categorical = FALSE,prob = .5,plot=FALSE))
pFc<-pFc[[1]]+ggthemes::theme_few()+scale_color_viridis_d()+scale_fill_viridis_d()+scale_x_continuous(breaks=c(0,1),labels = labs2)+ylab("")+xlab("Daytime temperature")+ylim(0,1)

pPc<-plot(conditional_effects(mod.likePc, "Light:phase", categorical = FALSE,prob = .5,plot=FALSE))
pPc<-pPc[[1]]+ggthemes::theme_few()+scale_color_viridis_d()+scale_fill_viridis_d()+scale_x_continuous(breaks=c(0,1),labels = labs3)+ylab("")+xlab("Photoperiod")+ylim(0,1)


ggpubr::ggarrange(ggpubr::ggarrange(p11,p22,p33,ncol=3,common.legend = TRUE),ggpubr::ggarrange(pCc,pFc,pPc,ncol=3,common.legend = TRUE),common.legend = TRUE,ncol=1)






dicho<-gather(dat,"phase","day",5:6)

comperDich<-filter(dicho,GEN.SPA=="COM.PER")
comperDich<-filter(comperDich,name!="COMPER ? HF DB") 

phen.comper.ind<-brm(day~phase*Force*Light+(phase*Force*Light|name),
                 warmup=3000,iter=4000,control=list(adapt_delta=.95),
                 data=comperDich) #40 obs

phen.comper<-brm(day~phase*Force*Light+(1|name),
                     warmup=3000,iter=4000,control=list(adapt_delta=.95),
                     data=comperDich) #40 obs

longday<-filter(comperDich,Light==1)
phen.comper2<-brm(day~phase*Force+(1|name),
                 warmup=3000,iter=4000,control=list(adapt_delta=.95),
                 data=longday) #40 obs

fixef(phen.comper2)
conditional_effects(phen.comper2,prob=.95)



nd.dich<-data.frame(Force=c(0,1,0,1),Light=c(0,0,1,1),phase=rep(c("flo_dayF","flo_dayM"),each=4))
namer<-data.frame(name=rep(unique(comperDich$name),each=8))
nd.dich<-cbind(namer,nd.dich)


plotter<-epred_draws(phen.comper.ind,newdata = nd.dich,ndraws = 100)
plotter$grouper<-paste(plotter$phase,plotter$.draw,plotter$name)
daynames<-as_labeller(c(`0`="8 hours",`1`="12 hours"))

ggplot(plotter,aes(Force,.epred))+stat_pointinterval(aes(color=phase),.width = c(.5,.89))+
       geom_line(aes(color=phase,group=grouper),size=0.02)+
  facet_grid(Light~name,labeller=daynames)+ggthemes::theme_few()+scale_color_viridis_d(labels=c("pistlate","staminate"))+
  ylab("Day of flowering")+scale_x_continuous(name="Daytime temperature",breaks=c(0,1),labels = c(18,24))+coord_cartesian(ylim=c(1,50))


nd.dich2<-data.frame(Force=c(0,1,0,1),Light=c(0,0,1,1),phase=rep(c("flo_dayF","flo_dayM"),each=4))
plotter2<-epred_draws(phen.comper,newdata = nd.dich2,ndraws = 1000,re_formula = NA)
plotter2$grouper<-paste(plotter2$phase,plotter2$.draw)


ggplot(plotter2,aes(Force,.epred))+
  geom_line(aes(color=phase,group=grouper),size=0.01)+stat_pointinterval(aes(color=phase),.width = c(.5,.89))+
  facet_wrap(~Light,labeller=daynames)+ggthemes::theme_few()+scale_color_viridis_d(labels=c("pistlate","staminate"))+
  ylab("Day of flowering")+scale_x_continuous(name="Daytime temperature",breaks=c(0,1),labels = c(18,24))+coord_cartesian(ylim=c(15,40),xlim=c(0,1))



nd.clim<-data.frame(Force=c(0,.66),Light=c(1,.875),phase=rep(c(scale_color_viridis_d()nd.clim<-data.frame(Force=c(0,.66),Light=c(1,.875),phase=rep(c("flo_dayF","flo_dayM"),each=2))
climpred<-fitted(phen.comper,newdata = nd.clim,probs = c(.25,.75),re_formula = NA)
climpred<-cbind(nd.clim,climpred)
climpred$scenario<-ifelse(climpred$Force==.66 &climpred$Light==.875,1,0)


ggplot(climpred,aes(scenario,Estimate))+
  geom_point(aes(color=phase),position=pd,size=4)+geom_line(aes(color=phase),linetype="dashed")+
  geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=phase),width=0,position=pd)+
  ggthemes::theme_few()+scale_color_viridis_d(labels=c("pistlate","staminate"))+
  ylab("Day of flowering")+scale_x_continuous(name="Climate scenario",breaks=c(0,1),labels = c("base","warming"))





    



otro<-dat %>% group_by(name,GEN.SPA,treatment,Force,Light,Chill) %>% dplyr::summarise(mean_m=mean(flo_dayM,na.rm=TRUE),mean_f=mean(flo_dayF,na.rm=TRUE))
otro<-filter(otro,GEN.SPA=="COM.PER")
otro$dichogamy<- otro$mean_m-otro$mean_f

phen.order<-brm(dichogamy~Force*Light+(1|name),iter=4000,warmup=3000,control=list(adapt_delta=0.95),data=otro)

conditional_effects(phen.order,prob=.75)















dicho2<datdicho2<-gather(dicho,"phase","day",7:8)


ggplot(dicho2,aes(as.factor(Light),day))+stat_summary(aes(color=phase))+facet_wrap(GEN.SPA~Force,scales="free")
ggplot(dicho,aes(as.factor(Force),dichogamy))+stat_summary()+facet_wrap(GEN.SPA~Light,scales="free")+geom_hline(yintercept=0)


mod1<-brm(dichogamy~Force*Light*GEN.SPA,data=dicho)
summary(mod1)

comper2<-filter(dicho2,GEN.SPA=="COM.PER")
comper<-filter(dicho,GEN.SPA=="COM.PER")

mod2<-brm(day~phase*Force*Light,data=comper2)

mod2a<-brm(dichogamy~Force*Light,data=comper)

nd.dich<-data.frame(Force=c(0,1,0,1),Light=c(0,0,1,1))
predo<-fitted(mod2a,newdata = nd.dich,probs = c(.25,.75))
predo<-cbind(nd.dich,predo)
ggplot(predo,aes(as.factor(Force),Estimate))+geom_point()+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`),width=.1)+facet_wrap(~Light,scales="free")+geom_hline(yintercept=0)




fixef(mod2,probs = c(.25,.75))
fixef(mod2a,probs = c(.25,.75))

nd<-data.frame(Force=c(0,1,0,1),Light=c(0,0,1,1),phase=rep(c("mean_m","mean_f"),each=4))
nd$GEN.SPA<-"COM.PER"
nd2<-nd

nd2$GEN.SPA<-"COR.COR"


pred<-fitted(mod2,newdata = nd,probs = c(.3,.7))
pred<-cbind(nd,pred)

ggplot(pred,aes(as.factor(Force),Estimate))+geom_point(aes(color=phase,shape=phase))+geom_errorbar(aes(ymin=`Q30`,ymax=`Q70`,color=phase),width=0)+facet_wrap(~Light,scales="free")
