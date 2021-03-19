###  extracts and add posteriors to group and plot posteriors by xls classes
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library("tidybayes")
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("~/Documents/git/proterant/FLOBUDS")
load("writing/flobud.main.mods.Rda")

goober<-mod.flo.int %>%
  spread_draws(b_Intercept,b_Force,b_Chill,b_Light,`b_Chill:Light`,`b_Chill:Force`,`b_Light:Force`)

goober2<-mod.flo.int%>%
  spread_draws(r_GEN.SPA[GEN.SPA,term])%>%
  spread(term,r_GEN.SPA) 
colnames(goober2)
flooby<-left_join(goober2,goober,by=c(".chain",".iteration",".draw")) %>% 
  mutate(Intercept_real=b_Intercept+Intercept,
         Force_real=+b_Force+Force,
         Chill_real=b_Chill+Chill,
         Light_real=+b_Light+Light,
         ChillxLight_real=`b_Chill:Light`+`Chill:Light`,
         ChillxForce_real=`b_Chill:Force`+`Chill:Force`,
         LightxForce_real=`b_Light:Force` + `Light:Force`)%>%
  mutate(basepred=Intercept_real+Light_real*1+Chill_real*.67+
           Force_real*0+ChillxLight_real*.67+ChillxForce_real*0+
           LightxForce_real*0) %>%
  mutate(warm_only=Intercept_real+Light_real*1+Chill_real*.67+
           Force_real*1+ChillxLight_real*.67+ChillxForce_real*.67+
           LightxForce_real*1)%>%
  mutate(more_chill=Intercept_real+Light_real*1+Chill_real*1+
           Force_real*1+ChillxLight_real*1+ChillxForce_real*1+
           LightxForce_real*1)%>%
  mutate(less_chill=Intercept_real+Light_real*1+Chill_real*0+
           Force_real*1+ChillxLight_real*0+ChillxForce_real*0+
           LightxForce_real*1)%>%
  gather(scenario,dof,26:29)


flooby$category<-NA
flooby$category[which(flooby$GEN.SPA %in% c("COR.COR","COM.PER","ACE.RUB"))]<-"flowering-first"
flooby$category[which(flooby$GEN.SPA %in% c("VAC.COR","ILE.MUC","ACE.PEN"))]<-"concurrent"
flooby$category[which(flooby$GEN.SPA %in% c("PRU.PEN","ILE.VER","PRU.VIR","VIB.ACE"))]<-"leafing-first"

flooby<-flooby%>% group_by(category,scenario,GEN.SPA)%>%
  mean_qi(dof,.width=0.5)


ggplot(flooby,aes(scenario,dof))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper))+facet_wrap(~category)

loober<-mod.lo.int %>%
  spread_draws(b_Intercept,b_Force,b_Chill,b_Light,`b_Chill:Light`,`b_Chill:Force`,`b_Light:Force`)

loober2<-mod.lo.int%>%
  spread_draws(r_GEN.SPA[GEN.SPA,term])%>%
  spread(term,r_GEN.SPA) 

looby<-left_join(loober2,loober,by=c(".chain",".iteration",".draw")) %>% 
  mutate(Intercept_real=b_Intercept+Intercept,
         Force_real=+b_Force+Force,
         Chill_real=b_Chill+Chill,
         Light_real=+b_Light+Light,
         ChillxLight_real=`b_Chill:Light`+`Chill:Light`,
         ChillxForce_real=`b_Chill:Force`+`Chill:Force`,
         LightxForce_real=`b_Light:Force` + `Light:Force`)%>%
  mutate(basepred=Intercept_real+Light_real*1+Chill_real*.67+
           Force_real*0+ChillxLight_real*.67+ChillxForce_real*0+
           LightxForce_real*0) %>%
  mutate(warm_only=Intercept_real+Light_real*1+Chill_real*.67+
           Force_real*1+ChillxLight_real*.67+ChillxForce_real*.67+
           LightxForce_real*1)%>%
  mutate(more_chill=Intercept_real+Light_real*1+Chill_real*1+
           Force_real*1+ChillxLight_real*1+ChillxForce_real*1+
           LightxForce_real*1)%>%
  mutate(less_chill=Intercept_real+Light_real*1+Chill_real*0+
           Force_real*1+ChillxLight_real*0+ChillxForce_real*0+
           LightxForce_real*1)%>%
  gather(scenario,dof,26:29)

looby$category<-NA
looby$category[which(looby$GEN.SPA %in% c("COR.COR","COM.PER","ACE.RUB"))]<-"flowering-first"
looby$category[which(looby$GEN.SPA %in% c("VAC.COR","ILE.MUC","ACE.PEN"))]<-"concurrent"
looby$category[which(looby$GEN.SPA %in% c("PRU.PEN","ILE.VER","PRU.VIR","VIB.ACE"))]<-"leafing-first"

looby<-looby%>% group_by(category,scenario,GEN.SPA)%>%
  mean_qi(dof,.width=0.5)


#######
bboober<-mod.bb.int %>%
  spread_draws(b_Intercept,b_Force,b_Chill,b_Light,`b_Chill:Light`,`b_Chill:Force`,`b_Light:Force`)

bboober2<-mod.bb.int%>%
  spread_draws(r_GEN.SPA[GEN.SPA,term])%>%
  spread(term,r_GEN.SPA) 

bbooby<-left_join(bboober2,bboober,by=c(".chain",".iteration",".draw")) %>% 
  mutate(Intercept_real=b_Intercept+Intercept,
         Force_real=+b_Force+Force,
         Chill_real=b_Chill+Chill,
         Light_real=+b_Light+Light,
         ChillxLight_real=`b_Chill:Light`+`Chill:Light`,
         ChillxForce_real=`b_Chill:Force`+`Chill:Force`,
         LightxForce_real=`b_Light:Force` + `Light:Force`)%>%
  mutate(basepred=Intercept_real+Light_real*1+Chill_real*.67+
           Force_real*0+ChillxLight_real*.67+ChillxForce_real*0+
           LightxForce_real*0) %>%
  mutate(warm_only=Intercept_real+Light_real*1+Chill_real*.67+
           Force_real*1+ChillxLight_real*.67+ChillxForce_real*.67+
           LightxForce_real*1)%>%
  mutate(more_chill=Intercept_real+Light_real*1+Chill_real*1+
           Force_real*1+ChillxLight_real*1+ChillxForce_real*1+
           LightxForce_real*1)%>%
  mutate(less_chill=Intercept_real+Light_real*1+Chill_real*0+
           Force_real*1+ChillxLight_real*0+ChillxForce_real*0+
           LightxForce_real*1)%>%
  gather(scenario,dof,26:29)

bbooby$category<-NA
bbooby$category[which(bbooby$GEN.SPA %in% c("COR.COR","COM.PER","ACE.RUB"))]<-"flowering-first"
bbooby$category[which(bbooby$GEN.SPA %in% c("VAC.COR","ILE.MUC","ACE.PEN"))]<-"concurrent"
bbooby$category[which(bbooby$GEN.SPA %in% c("PRU.PEN","ILE.VER","PRU.VIR","VIB.ACE"))]<-"leafing-first"

bbooby<-bbooby%>% group_by(category,scenario,GEN.SPA)%>%
  mean_qi(dof,.width=0.5)
#ggplot(flooby,aes(scenario,dof))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper))+facet_wrap(~category)


looby$phase<-"leafout"
bbooby$phase<-"budburst"
flooby$phase="flowering"

posty<-rbind(flooby,looby,bbooby)
posty$class<-paste(posty$phase,posty$category)
posty$pollination<-ifelse(posty$category=="flowering-first","abotic","biotic")

posty$scenario2 = factor(posty$scenario, levels=c('basepred',"warm_only","less_chill","more_chill"))

posty$scenario3<-NA

posty$scenario3[which(posty$scenario2=="basepred")]<-"base"
posty$scenario3[which(posty$scenario2=="warm_only")]<-"forcing"
posty$scenario3[which(posty$scenario2=="less_chill")]<-"updown"
posty$scenario3[which(posty$scenario2=="more_chill")]<-"upup"




#setEPS()
#postscript("Plots/postergroups2.eps",width = 8, height = 4)

pd=position_dodge(width=1.5)
one<-posty %>%
  arrange(dof) %>%
  mutate(category = factor(category, levels=c("flowering-first","concurrent","leafing-first"))) %>%
  ggplot(aes(scenario3,dof,shape=phase,color=phase))+
  geom_errorbar(aes(ymin=.lower,ymax=.upper),width=0.0)+
  geom_point(size=2.5)+
  ggthemes::theme_few()+
  facet_wrap( ~ GEN.SPA,nrow=1,switch = "x")+
  theme(strip.background = element_blank(),
  strip.text.x = element_blank())+
  scale_color_manual(values = c("lightgreen","purple","darkgreen"))+
  scale_shape_manual(values=c(17,19,15))+
  scale_x_discrete(labels=c("Baseline",
                            "+Forcing",
                            "+Force/-Chill",
                            "+Force/+Chill"))+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_text(angle=270,vjust=0),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        strip.placement = "bottom",
        panel.spacing.x=unit(0,"cm"),
        panel.spacing.y=unit(0,"cm"))+
  ylim(0,110)+
  ylab("Day of phenological event")+
  theme(legend.position = "bottom")
 
  


one
dev.off()  

posto<-spread(posty,phase,dof)

posto2<-posto%>%group_by(GEN.SPA,scenario3) %>%summarise(interphase=mean(leafout,na.rm=TRUE)-mean(flowering,na.rm=TRUE))
posto3<-posto%>%group_by(GEN.SPA,scenario3) %>%summarise(interphase=mean(budburst,na.rm=TRUE)-mean(flowering,na.rm=TRUE))

loerror<-filter(posto,!is.na(leafout))

ferror<-filter(posto,!is.na(flowering))
bberror<-filter(posto,!is.na(budburst))



posto2$interphase2<-abs(posto2$interphase)
posto3$interphase2<-abs(posto3$interphase)

posto2$Hyst<-ifelse(posto2$interphase>=0,"Flower-first", "Vegetative-First ")
posto3$Hyst<-ifelse(posto3$interphase>=0,"Flower-first", "Vegetative-First ")

posto2$CI<-abs(loerror$.lower-ferror$.upper)*.5
posto3$CI<-abs(bberror$.lower-ferror$.upper)*.5

two<-posto2 %>%
  arrange(interphase) %>%
  ggplot()+geom_point(aes(scenario3,y=interphase),position=pd,size=2,shape=8)+
  geom_errorbar(aes(scenario3,ymin=interphase-CI,ymax=interphase+CI),width=0)+
  #geom_text(aes(scenario3,y=interphase2+CI+3,label=paste(round(interphase2,digits=0),"(",round(CI,digits=0),")"),color=Hyst),position=pd,size=3)+
  ggthemes::theme_few(base_size = 10)+facet_wrap( ~ GEN.SPA,nrow=1,scales="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+
  ylab("Flowering-leafout \ninterphase (days)")+
 scale_color_manual(values=c("black","red"))+
  scale_x_discrete(labels=c("Baseline",
                            "+Forcing",
                            "+Force/-Chill",
                            "+Force/+Chill"))+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+theme(legend.position ="none")+
  geom_hline(yintercept=0,linetype="dotted")
 

two

posto3$species<-NA
posto3$species[which(posto3$GEN.SPA=="ACE.PEN")]<-"Acer pensylvanicum"
posto3$species[which(posto3$GEN.SPA=="ACE.RUB")]<-"Acer rubrum"
posto3$species[which(posto3$GEN.SPA=="COM.PER")]<-"Comptonia peregrina"
posto3$species[which(posto3$GEN.SPA=="COR.COR")]<-"Corylus cornuta"
posto3$species[which(posto3$GEN.SPA=="ILE.MUC")]<-"Ilex mucronata"
posto3$species[which(posto3$GEN.SPA=="ILE.VER")]<-"Ilex verticillata"
posto3$species[which(posto3$GEN.SPA=="PRU.PEN")]<-"Prunus pensylvanica"
posto3$species[which(posto3$GEN.SPA=="PRU.VIR")]<-"Prunus virginiana"
posto3$species[which(posto3$GEN.SPA=="VAC.COR")]<-"Vaccinium corymbosum"
posto3$species[which(posto3$GEN.SPA=="VIB.ACE")]<-"Viburnum acerifolium"

  three<-posto3 %>%
    arrange(interphase) %>%
    ggplot()+geom_point(aes(scenario3,y=interphase),position=pd,size=2,shape=8)+
    geom_errorbar(aes(scenario3,ymin=interphase-CI,ymax=interphase+CI),width=0)+
    #geom_text(aes(scenario3,y=interphase2+4,label=paste(round(interphase2,digits=0),round(CI,digits=2)),color=Hyst),position=pd,size=3)+
    ggthemes::theme_few(base_size = 10)+facet_wrap( ~ species,nrow=1,scales="free_y")+
   
    ylab("Flowering-budburst \ninterphase (days)")+
   scale_color_manual(values=c("black","red"))+
    scale_x_discrete(labels=c("Baseline",
                              "+Forcing",
                              "+Force/-Chill",
                              "+Force/+Chill"))+
    theme(axis.title.x =element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank())+theme(legend.position ="none")+
     theme(strip.text = element_text(face = "italic"))+geom_hline(yintercept=0,linetype="dotted")

  png("Plots/species_projections.png",width = 14, height= 10,units = 'in',res = 300)  
  ggpubr::ggarrange(three,two,one,ncol=1,nrow=3,heights = c(.2,.2,1))
dev.off()  
  

