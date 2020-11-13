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

loober<-mod.bb.int %>%
  spread_draws(b_Intercept,b_Force,b_Chill,b_Light,`b_Chill:Light`,`b_Chill:Force`,`b_Light:Force`)

loober2<-mod.bb.int%>%
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

ggplot(flooby,aes(scenario,dof))+geom_point()+geom_errorbar(aes(ymin=.lower,ymax=.upper))+facet_wrap(~category)
looby$phase<-"leaf"
flooby$phase="flower"

posty<-rbind(flooby,looby)
posty$class<-paste(posty$phase,posty$category)
posty$pollination<-ifelse(posty$category=="flowering-first","abotic","biotic")

posty$scenario2 = factor(posty$scenario, levels=c('basepred',"warm_only","less_chill","more_chill"))

posty$scenario3<-NA

posty$scenario3[which(posty$scenario2=="basepred")]<-"base prediction"
posty$scenario3[which(posty$scenario2=="warm_only")]<-"warming only"
posty$scenario3[which(posty$scenario2=="less_chill")]<-"decreasing chilling"
posty$scenario3[which(posty$scenario2=="more_chill")]<-"increasing chilling"
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
  scale_shape_manual(values=c(19,17))+
  scale_x_discrete(labels =c("FL","F=L","LF"))
dev.off()  

#setEPS()
#postscript("Plots/postergroups2.eps",width = 8, height = 4)
png("Plots/posteriorgroups_go.png",width = 7, height= 2.5,units = 'in',res = 200)
pd=position_dodge(width=1.5)
posty %>%
  arrange(dof) %>%
  mutate(category = factor(category, levels=c("flowering-first","concurrent","leafing-first"))) %>%
  ggplot(aes(0,dof,color=category,shape=phase,group =category))+
  geom_errorbar(aes(ymin=.lower,ymax=.upper),position=pd,width=0.4,alpha=0.5)+
  geom_point(size=2.5,position=pd)+
  ggthemes::theme_few(base_size = 11)+
  facet_wrap( ~ scenario3,drop=TRUE, scales="free_x",switch="x",ncol=4)+
  theme(axis.title.x =element_blank(),
    strip.placement = "bottom",
    strip.background = element_rect(color="black"),
    panel.spacing.x=unit(0,"cm"),
    panel.spacing.y=unit(0,"cm"))+
  ylim(5,100)+
  ylab("Day of phenological event")+
  scale_color_manual(values = c("#fc8d62","#8da0cb","#66c2a5"))+
  scale_shape_manual(values=c(19,17))+
  scale_x_discrete(labels =c("FL","F=L","LF"))
dev.off()  
