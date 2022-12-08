rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(gtable)
library(ggthemes)

library(tidyverse)





x<-c(320,400,440,520)
y<-c(8,8,16,16)
 
photoperiod<-c("low","low","high","high")
z <- (y*-.2)+(x*-.04)+50
forcing<-c(0,1,0,1)
dat<-data.frame(x,y,z,photoperiod,forcing)

ploa<-ggplot(dat,aes(x,z))+stat_smooth(method="lm",fullrange = FALSE,aes(color=photoperiod,linetype="estimated effect"),size=1.5)+stat_smooth(method="lm",fullrange = TRUE,size=0.5,aes(color=photoperiod,linetype="true effect"))+scale_color_viridis_d(option="turbo")+ggthemes::theme_few()+ylab("Day of budburst")+
  xlab("Thermal sums")+#+ylim(34,48)+
scale_linetype_manual(name = "effect",
                   values = c( "estimated effect" = "solid", "true effect" = "dotted"))+

 annotate("text", x = 470, y = 36, label = "no interaction")


plob<-ggplot(dat,aes(forcing,z))+stat_smooth(method="lm",fullrange = FALSE,aes(color=photoperiod),size=1.5)+scale_color_viridis_d(option="turbo")+ggthemes::theme_few()+ylab("Day of budburst")+
  annotate("text", x = .9, y = 36, label = "no interaction")+
  scale_x_continuous(limits = c(-.25,1.25),breaks=0:1,labels=c(labels=c(expression(""*20/10~degree*C)),expression(""*25/15~degree*C)))+xlab("Forcing treatment")#+
  
ggpubr::ggarrange(ploa,plob)

z2 <- y*-.2+x*-.04+(x*y*.001)+80


ploc<-ggplot(dat,aes(x,z2))+stat_smooth(method="lm",aes(color=photoperiod),size=1.5)+stat_smooth(method="lm",fullrange = TRUE,linetype="dotted",size=0.5,aes(color=photoperiod))+
  scale_color_viridis_d(option="turbo")+ggthemes::theme_few()+ylab("Day of budburst")+
  xlab("Thermal sums")+ annotate("text", x = 470, y = 70, label = "with interaction")

plod<-ggplot(dat,aes(forcing,z2))+stat_smooth(method="lm",fullrange = FALSE,aes(color=photoperiod),size=1.5)+scale_color_viridis_d(option="turbo")+ggthemes::theme_few()+ylab("Day of budburst")+
  scale_x_continuous(limits = c(-.25,1.25),breaks=0:1,labels=c(c(expression(""*20/10~degree*C)),expression(""*25/15~degree*C)))+xlab("Forcing treatment")+
  annotate("text", x = .9, y = 69, label = "with interaction")

?jpeg()
jpeg("~/git/proterant/FLOBUDS/Plots/periodicity_figures/apparent.jpeg",width = 8,height=6,unit="in",res=200)
ggpubr::ggarrange(ploa,plob,ploc,plod,common.legend = TRUE,legend = "right",labels  = c("a)","b)","c)","d)"))
dev.off()


(25*8)+(15*16)
(20*8)+(10*16)

force<-c(360,480,360,480,320,400,440,520)
photo<-c(8,8,16,16,8,8,16,16)
photoperiod<-c("low","low","high","high","low","low","high","high")


dat2<-data.frame(force,photo,photoperiod)

dat2$budburst<-(force*-.07)+(photo*-.04)+100
dat2$forcing<-c(0,1,0,1,0,1,0,1)

dat2$design<-c("orthoginal","orthoginal","orthoginal","orthoginal","exp. covariation","exp. covariation","exp. covariation","exp. covariation")

ggplot(dat2,aes(force,budburst,color=photoperiod,linetype=design))+geom_smooth(method="lm")+
  ggthemes::theme_few()+scale_color_viridis_d(begin = .1,end = .5,direction = -1)+
  ylab("day of budburst")+scale_linetype_manual(values =c("dotdash","solid"))+facet_wrap(~design)

ploty1<-ggplot(dat2,aes(forcing,budburst,color=photoperiod,linetype=design))+geom_smooth(method="lm")+
  ggthemes::theme_few()+scale_x_continuous(breaks=c(0,1),labels = c("low","high"))+scale_color_viridis_d(begin = .1,end = .5,direction = -1)+
  ylab("day of budburst")+scale_linetype_manual(values =c("dotdash","solid"))

dat2$budburst2<-(force*-.07)+(photo*-.04)+(force*photo*-.0007)+100
dat2$budburst3<-(force*-.07)+(photo*-.04)+(force*photo*.0007)+100
ploty2<-ggplot(dat2,aes(forcing,budburst2,color=photoperiod,linetype=design))+geom_smooth(method="lm")+
  ggthemes::theme_few()+scale_x_continuous(breaks=c(0,1),labels = c("low","high"))+
  scale_color_viridis_d(begin = .1,end = .5,direction = -1)+ylab("day of budburst")+
  scale_linetype_manual(values =c("dotdash","solid"))+theme(legend.position = "none")

ploty3<-ggplot(dat2,aes(forcing,budburst3,color=photoperiod,linetype=design))+geom_smooth(method="lm")+
  ggthemes::theme_few()+scale_x_continuous(breaks=c(0,1),labels = c("low","high"))+
  scale_color_viridis_d(begin = .1,end = .5,direction = -1)+ylab("day of budburst")+
  scale_linetype_manual(values =c("dotdash","solid"))+theme(legend.position = "none")

inters<-ggpubr::ggarrange(ploty2,ploty3,labels=c("b)","c)"))
ggpubr::ggarrange(ploty1,inters,ncol=1,nrow=2,common.legend=TRUE,labels=c("a)"))

dat1<-filter

fig <- plot_ly(x = ~x, y = ~y, z = ~z,color=~rev(cats), colors = c( '#0C4B8E','#BF382A'), type = 'mesh3d',size=0.001)


fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Photoperiod'),
                                   yaxis = list(title = 'Forcing'),
                                   zaxis = list(title = 'Day of Budburst')))

fig










library(ggplot2)
Temperature<- c(25,15,25,15)
Photoperiod<-c(8,8,16,16)
Treatments<-c("short/high","short/low","long/high","long/low")
ortho<-data.frame(Temperature,Photoperiod,Treatments)

one<-ggplot(ortho,aes(Photoperiod,Temperature))+geom_rect(xmin=8,xmax=16,ymin=15,ymax=25,alpha=0.1,fill="royalblue2")+
  geom_point(aes(shape=Treatments),size=3)+ylim(10,30)+xlim(6,18)+
  scale_shape_manual(values=c(0,1,2,5))+ggthemes::theme_few()

?scale_shape_manual()
Temperature<- c(25,13,28,15)
Photoperiod<-c(8,8,16,16)
Treatments<-c("short/high","short/low","long/high","long/low")

noortho<-data.frame(Temperature,Photoperiod,Treatments)
two<-ggplot(noortho,aes(Photoperiod,Temperature))+geom_rect(xmin=8,xmax=16,ymin=15,ymax=25,alpha=0.1,fill="royalblue2")+
  geom_point(aes(shape=Treatments),size=3)+
  ylim(10,30)+xlim(6,18)+scale_shape_manual(values=c(0,1,2,5))+ggthemes::theme_few()


Temperature<- c(25,15,25)
Photoperiod<-c(8,8,16)
Treatments<-c("short/high","short/low","long/high")
noint<-data.frame(Temperature,Photoperiod,Treatments)
three<-ggplot(noint,aes(Photoperiod,Temperature))+geom_rect(xmin=8,xmax=16,ymin=15,ymax=25,alpha=0.1,fill="royalblue2")+
  geom_point(aes(shape=Treatments),size=3)+
  ylim(10,30)+xlim(6,18)+scale_shape_manual(values=c(0,1,2))+ggthemes::theme_few()

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/factorial.jpeg",width = 6,height=8,unit="in",res=300)
ggpubr::ggarrange(one,two,three,labels=c("a)","b)","c)"),ncol=1,nrow=3,common.legend = TRUE,legend= "right")
dev.off()

GDH<-c(7.5,12.5,20,22.5,10,20,20,10)

photo<-c(11,11,15,15,11,11,15,15)
design<-c("non-orthoginal","non-orthoginal","non-orthoginal","non-orthoginal","orthoginal","orthoginal","orthoginal","orthoginal") 

dat<-data.frame(GDH,photo,design)
dat$phen<--1.9*dat$GDH-1.2*dat$photo +1000

x<-c(10.8,10.8,15.2,15.2) 
y<-c(382,482,382,482)

dat.jr<-filter(dat,design=="non-orthoginal")

ggplot(dat.jr,aes(x=photo,y=GDH,z=phen))+theme_void()+axes_3D()+
  stat_wireframe(alpha=.5,phi=30)

xy<-data.frame(x,y)





ploty<-ggplot(dat,(aes(x=photo)))+geom_point(aes(y=GDH,color=design),size=4)+xlim(8,18)+ylim(350,520)+theme_bw()+geom_polygon(aes(x=photo,y=GDH,fill=design),alpha=0.5)+scale_fill_manual(values=c("tomato1","royalblue2"))

plotty2<-ploty+geom_rect(dat,mapping=aes(xmin=10.9,xmax=15.1,ymin=382,ymax=482),fill=NA,color="black",size=0.2,linetype="dashed")

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/orthog.jpeg",width = 11,height=8,unit="in",res=300)
plotty2+geom_point(xy,mapping=aes(x,y),shape=8,size=4)+theme(legend.position = "left")+ylab("Temperature (thermal sums")+xlab("Photoperiod")
dev.off()


###try 3d


ggplot(dat,aes(x=photo,y=GDH,z=phen,color=design,,fill=design))+theme_void()+axes_3D()+
  stat_wireframe(alpha=.5,phi=30)
 #+stat_3D(geom="path")



photoperiod<-rep(seq(8,17, by=0.1),10)
forcing<-rep(seq(0,9,by=0.1),each=10)


dat2<-data.frame(photoperiod=photoperiod,forcing=forcing)






  
  
dat2$phen<--2.5*dat2$forcing-1.5*dat2$photoperiod+100
dat3<-dat2

dat3$forcing<-rep(c(seq(-2,7,by=0.1),seq(2,11,by=0.1)),each=5)

dat3$phen<--2.5*dat3$forcing-1.5*dat3$photoperiod+0.2+100
dat2$design<-"non-covarying"
dat3$design<-"covarying"
dat4<-rbind(dat3,dat2)

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/orthog.jpeg",width = 8,height=8,unit="in",res=200)
ggplot(dat4, aes(photoperiod,forcing,z=phen)) +
  axes_3D() + theme_void()+
  labs_3D(labs=c("photoperiod", "forcing", "Phenological event"))+scale_color_viridis_d()+

stat_wireframe(alpha=0.6,aes(color=design))
dev.off()


ggplot(dat3, aes(photoperiod,forcing,z=phen)) +theme_void()+
  axes_3D(phi=30) +  stat_wireframe(alpha=.5,phi=30)


  theme_void()

?(gg3d)
########
temp<-c(12,12,12,12)
light<-c(11,11,15,15)
treat<-c("warm/short","cool/short","warm/long","warm/long")
dawn<-data.frame(temp,light,treat)
dawn$dawn.offset<-dawn$temp-dawn$light



xx<-seq(0,24,by=0.1)
yy<-rep(20,241)
zz<-rep(15,241)
aa<-c(rep(c(15,15,15,15,15,15,15,15),each=10),rep(c(25,25,25,25,25,25,25,25),each=10),rep(c(15,15,15,15,15,15,15,15),each=10),15)
bb<-c(rep(c(15,15,15,15,15,15),each=10),rep(c(25,25,25,25,25,25,25,25,25,25,25,25),each=10),rep(c(15,15,15,15,15,15),each=10),15)
cc<-c(rep(c(10,10,10,10,10,10,10,10),each=10),rep(c(20,20,20,20,20,20,20,20),each=10),rep(c(10,10,10,10,10,10,10,10),each=10),10)
dd<-c(rep(c(10,10,10,10,10,10),each=10),rep(c(20,20,20,20,20,20,20,20,20,20,20,20),each=10),rep(c(10,10,10,10,10,10),each=10),10)

ee<-c(rep(c(15,15,15,15,15,15,15),each=10),rep(c(25,25,25,25,25,25,25,25,25,25,25,25),each=10),rep(c(15,15,15,15,15),each=10),15)
ff<-c(rep(c(10,10,10,10,10,10,10),each=10),rep(c(20,20,20,20,20,20,20,20,20,20,20,20),each=10),rep(c(10,10,10,10,10),each=10),10)


gg<-c(22,22,22,22,22,22,22,22,31,31,31,31,31,31,31,22,22,22,22,22,22,22,22,22)
hh<-c(24,24,24,24,24,24,26,26,26,26,26,26,26,26,26,26,26,26,24,24,24,24,24,24)
ii<-c(17,17,17,17,17,17,17,17,26,26,26,26,26,26,26,26,17,17,17,17,17,17,17,17)
jj<-c(19,19,19,19,19,19,21,21,21,21,21,21,21,21,21,21,21,21,19,19,19,19,19,19)

dat.simple<-data.frame(xx,yy,zz,aa,bb,cc,dd,ee,ff)#,gg,hh,ii)
a<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 20 \n, diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
b<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=4,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T= 20 \n diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
c<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short\n daily mean T= 15 \n diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
d<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=4,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T= 15 \n diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
flat<-ggpubr::ggarrange(a, b,c, d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

aaa<-ggplot(dat.simple,aes(xx,aa))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,aa),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 23.3 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=4,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,bb),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T = 25 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,cc),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short \n daily mean T= 13.3 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=4,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,dd),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T=15 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))




jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/basic.jpeg",width = 8,height=8,unit="in",res=200)
ggarrange(ccc, ddd, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()



aaaa<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 20 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=4,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T= 20 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
cccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short \n daily mean T= 15 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
dddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=4,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T= 15 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
noncovarying<-ggpubr::ggarrange(aaaa, bbbb,cccc, dddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 


fully<-ggpubr::ggarrange(aaaa,a, bbbb,b,cccc, c,dddd,d, ncol=4, nrow=2, common.legend = TRUE, legend="bottom") 

aaaaa<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,gg),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 25 \n diurnal differnce= 9")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbbbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,hh),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T= 25 \n diurnal differnce= 2")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ccccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,ii),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short \n daily mean T= 15 \n diurnal differnce= 9")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ddddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,jj),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T= 15 \n diurnal differnce= 2")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
noncovarying2<-ggarrange(aaaaa, bbbbb,ccccc, ddddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 

chico<-ggpubr::ggarrange(flat,noncovarying,labels=c("a)","b)"))

plot.list <- lapply(list(chico,fully), 
                    function(p) p + theme(plot.background = element_rect(color = "black")))

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/designs.jpeg",width = 11,height=11,unit="in",res=200)

ggpubr::ggarrange(plotlist = plot.list,ncol=1,nrow=2,labels = c("","c)"))

grid.rect(width = 1, height = 0,.5, gp = gpar(lwd = 2, col = "black", fill = NA,hjust="left",vjust="topleft"))

grid.rect(width = 0,.5, height = 1, gp = gpar(lwd = 2, col = "black", fill = NA),y =1 )
annotate_figure()
dev.off()

