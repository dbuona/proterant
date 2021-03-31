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

flooby<-flooby%>% group_by(category,scenario)%>%
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

looby<-looby%>% group_by(category,scenario)%>%
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

bbooby<-bbooby%>% group_by(category,scenario)%>%
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



setEPS()
postscript("Plots/postergroups.eps",width = 8, height = 4)

pd=position_dodge(width=0.1)
posty %>%
  arrange(dof) %>%
  mutate(category = factor(category, levels=c("hyst","syn","ser"))) %>%
  ggplot(aes(category,dof,color=phase,shape=phase,fill=phase))+geom_point(size=2.5,position=pd)+geom_errorbar(aes(ymin=.lower,ymax=.upper),position=pd,width=0)+
  ggthemes::theme_few(base_size = 10)+
  facet_grid(. ~ scenario2, scales="free_x",switch="x",drop = TRUE)+
  theme(
    #axis.ticks.x=element_blank(),
    strip.placement = "bottom",
    strip.background = element_rect(color="black"),
    panel.spacing.x=unit(0,"cm"))+
  ylim(5,100)+
  ylab("Day of phenological event")+
  scale_color_manual(values=c("darkorchid3","darkgreen"))+
  scale_shape_manual(values=c(15,19,17))+
  scale_x_discrete(labels =c("FL","F=L","LF"))
dev.off()  

#setEPS()
#postscript("Plots/postergroups2.eps",width = 8, height = 4)

pd=position_dodge(width=1.5)
goo<-posty %>%
  arrange(dof) %>%
  mutate(category = factor(category, levels=c("flowering-first","concurrent","leafing-first"))) %>%
  ggplot(aes(0,dof,color=category,shape=phase,group =category))+
  geom_errorbar(aes(ymin=.lower,ymax=.upper),position=pd,width=0,alpha=0.5)+
  geom_point(size=2.5,position=pd)+
  ggthemes::theme_few(base_size = 11)+
  facet_wrap( ~ scenario3,drop=TRUE, scales="free_x",switch="x",ncol=4,
              labeller = as_labeller(c(`base` = "Baseline \n(C=6w,F=21°C)",`forcing` = "+Forcing \n(C=6w,F=+6°C)", `updown` = "+Force/-Chill \n(C=4w,F=+6°C)",`upup` = "+Force/+Chill \n(C=8w,F=+6°C)")))+
  theme(axis.title.x =element_blank(),
    strip.placement = "top",
    strip.background = element_rect(color="black"),
    panel.spacing.x=unit(0,"cm"),
    panel.spacing.y=unit(0,"cm"))+
  ylim(5,100)+
  ylab("Day of phenological event")+
  scale_color_manual(values = c("#fc8d62","#8da0cb","#66c2a5"))+
  scale_shape_manual(values=c(17,19,15))+
  scale_x_discrete(labels =c("FL","F=L","LF"))+
  theme(legend.position = "bottom")
  


dev.off()  

posto<-spread(posty,phase,dof)


posto2<-posto%>%group_by(category,scenario3) %>%summarise(interphase=mean(leafout,na.rm=TRUE)-mean(flowering,na.rm=TRUE))
posto3<-posto%>%group_by(category,scenario3) %>%summarise(interphase=mean(budburst,na.rm=TRUE)-mean(flowering,na.rm=TRUE))

loerror<-filter(posto,!is.na(leafout))

ferror<-filter(posto,!is.na(flowering))
bberror<-filter(posto,!is.na(budburst))
posto2$CI<-abs(loerror$.lower-ferror$.upper)*.5
posto3$CI<-abs(bberror$.lower-ferror$.upper)*.5

goo2<-posto2 %>%
  arrange(interphase) %>%
  mutate(category = factor(category, levels=c("flowering-first","concurrent","leafing-first"))) %>%
  ggplot()+geom_point(aes(0,y=interphase,color=category),position=pd,size=2,shape=8)+
  geom_errorbar(aes(0,ymin=interphase-CI,ymax=interphase+CI,color=category),position=pd,width=0)+
  #geom_text(aes(0,y=interphase+2,color=category,label=round(interphase,digits=1)),position=pd,size=3)+
  ggthemes::theme_few(base_size = 10)+facet_wrap( ~ scenario3,drop=TRUE,switch="x",ncol=4)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  theme(axis.title.x =element_blank(),
        strip.placement = "bottom",
        strip.background = element_rect(color="black"),
        panel.spacing.x=unit(0,"cm"),
        panel.spacing.y=unit(0,"cm"))+
  ylab("Flowering-leafout \ninterphase (days)")+
  
  scale_color_manual(values = c("#fc8d62","#8da0cb","#66c2a5"))+
  scale_shape_manual(values=c(17,19,15))+
  scale_x_discrete(labels =c("FL","F=L","LF"))+
  theme(legend.position ="none")+theme(legend.position = "none")+
  geom_hline(yintercept=0,linetype="dotted")
goo2  


goo3<-posto3 %>%
  arrange(interphase) %>%
  mutate(category = factor(category, levels=c("flowering-first","concurrent","leafing-first"))) %>%
  ggplot()+geom_point(aes(0,y=interphase,color=category),shape=8,position=pd,size=2)+
  geom_errorbar(aes(0,ymin=interphase-CI,ymax=interphase+CI,color=category),position=pd,width=0)+
  #geom_text(aes(0,y=interphase+2,color=category,label=round(interphase,digits=1)),position=pd,size=3)+
  ggthemes::theme_few(base_size = 10)+facet_wrap( ~ scenario3,drop=TRUE,switch="x",ncol=4)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  theme(axis.title.x =element_blank(),
        strip.placement = "bottom",
        strip.background = element_rect(color="black"),
        panel.spacing.x=unit(0,"cm"),
        panel.spacing.y=unit(0,"cm"))+
  ylab("Flowering-budburst \ninterphase (days)")+
  scale_color_manual(values = c("#fc8d62","#8da0cb","#66c2a5"))+
  scale_shape_manual(values=c(17,19,15))+
  scale_x_discrete(labels =c("FL","F=L","LF"))+
  theme(legend.position ="none")+theme(legend.position = "none")+
  geom_hline(yintercept=0,linetype="dotted")
goo3 


png("Plots/posteriorgroups_go.png",width = 8, height= 9,units = 'in',res = 200)
ggpubr::ggarrange(goo2,goo3,goo,nrow=3,ncol=1,heights = c(.2,.2,1),labels = c("a)","b)","c)"),label.x = .07)
dev.off()

### now do the
