x<-c(320,360,400,440,480,520)
y<-c(8,8,12,12,16,16)

photoperiod<-c("low","low","control","control","high","high")
z <- (y*-.2)+(x*-.04)+50
forcing<-c(0,1,0,1,0,1)
dat<-data.frame(x,y,z,photoperiod,forcing)

ploa<-ggplot(dat,aes(x,z))+stat_smooth(method="lm",fullrange = FALSE,aes(color=photoperiod,linetype="estimated effect"),size=1.5)+stat_smooth(method="lm",fullrange = TRUE,size=0.5,aes(color=photoperiod,linetype="true effect"))+scale_color_viridis_d(option="turbo")+ggthemes::theme_few()+ylab("Day of budburst")+
  xlab("Thermal sums")+#+ylim(34,48)+
  scale_linetype_manual(name = "effect",
                        values = c( "estimated effect" = "solid", "true effect" = "dotted"))

plob<-ggplot(dat,aes(forcing,z))+stat_smooth(method="lm",fullrange = FALSE,aes(color=photoperiod),size=1.5)+scale_color_viridis_d(option="turbo")+ggthemes::theme_few()+ylab("Day of budburst")+
  scale_x_continuous(limits = c(-.25,1.25),breaks=c(0,.5,1),labels=c(labels=c(expression(""*20/10~degree*C)),expression(""*22.5/12.5~degree*C),expression(""*25/15~degree*C)))+xlab("Forcing treatment")#+

ggpubr::ggarrange(ploa,plob,common.legend=TRUE,labels=c("a)","b)"))

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/apparent4spp.jpeg",width = 6,height=5,unit="in",res=200)
ggpubr::ggarrange(ploa,plob,common.legend=TRUE,labels=c("a)","b)"))
dev.off()


ggpubr::ggarrange(ploa,plob,common.legend = TRUE)
