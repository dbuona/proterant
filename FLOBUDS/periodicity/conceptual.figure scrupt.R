rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(gtable)
library(ggthemes)
library(gg3D)
library(tidyverse)



##FIg1
Temperature<- c(30,20,30,20)
Photoperiod<-c(8,8,12,12)
Treatments<-c("short/high","short/low","long/high","long/low")
ortho<-data.frame(Temperature,Photoperiod,Treatments)

one<-ggplot(ortho,aes(Photoperiod,Temperature))+geom_rect(xmin=8,xmax=12,ymin=20,ymax=30,alpha=0.1,fill="royalblue2")+
  geom_point(aes(shape=Treatments),size=3)+ylim(15,35)+xlim(6,14)+
  scale_shape_manual(values=c(0,1,2,5))+theme_few()

?scale_shape_manual()
Temperature<- c(30,17,33,20)
Photoperiod<-c(8,8,12,12)


noortho<-data.frame(Temperature,Photoperiod,Treatments)
two<-ggplot(noortho,aes(Photoperiod,Temperature))+geom_rect(xmin=8,xmax=12,ymin=20,ymax=30,alpha=0.1,fill="royalblue2")+
  geom_point(aes(shape=Treatments),size=3)+
  ylim(15,35)+xlim(6,14)+scale_shape_manual(values=c(0,1,2,5))+theme_few()


Temperature<- c(30,20,30)
Photoperiod<-c(8,8,12)
Treatments<-c("short/high","short/low","long/high")
noint<-data.frame(Temperature,Photoperiod,Treatments)
three<-ggplot(noint,aes(Photoperiod,Temperature))+geom_rect(xmin=8,xmax=12,ymin=20,ymax=30,alpha=0.1,fill="royalblue2")+
  geom_point(aes(shape=Treatments),size=3)+
  ylim(15,35)+xlim(6,14)+scale_shape_manual(values=c(0,1,2))+theme_few()

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/factorial.jpeg",width = 6,height=8,unit="in",res=300)
ggarrange(one,two,three,labels=c("a)","b)","c)"),ncol=1,nrow=3,common.legend = TRUE,legend= "right")
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



xx<-c(1:24)
yy<-rep(25,24)
zz<-rep(15,24)
aa<-c(20,20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,20,20,20,20,20,20,20,20)
bb<-c(20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20,20)
cc<-c(10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,10)
dd<-c(10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10)

ee<-c(20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20)
ff<-c(10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10)


gg<-c(22,22,22,22,22,22,22,22,31,31,31,31,31,31,31,22,22,22,22,22,22,22,22,22)
hh<-c(24,24,24,24,24,24,26,26,26,26,26,26,26,26,26,26,26,26,24,24,24,24,24,24)
ii<-c(17,17,17,17,17,17,17,17,26,26,26,26,26,26,26,26,17,17,17,17,17,17,17,17)
jj<-c(19,19,19,19,19,19,21,21,21,21,21,21,21,21,21,21,21,21,19,19,19,19,19,19)

dat.simple<-data.frame(xx,yy,zz,aa,bb,cc,dd,ee,ff,gg,hh,ii)
a<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 25 \n, diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
b<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T= 25 \n diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
c<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short\n daily mean T= 15 \n diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
d<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T= 15 \n diurnal differnce= NA")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
flat<-ggarrange(a, b,c, d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

aaa<-ggplot(dat.simple,aes(xx,aa))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,aa),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 23.3 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,bb),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T = 25 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,cc),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short \n daily mean T= 13.3 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,dd),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T=15 \n diurnal differnce= 10 C")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
covarying<-ggarrange(aaa, bbb,ccc, ddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 

aaaa<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 25 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T= 25 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
cccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short \n daily mean T= 15 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
dddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T= 15 \n diurnal differnce= 10")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
noncovarying<-ggarrange(aaaa, bbbb,cccc, dddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 


aaaaa<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,gg),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short \n daily mean T= 25 \n diurnal differnce= 9")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbbbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,hh),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long \n daily mean T= 25 \n diurnal differnce= 2")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ccccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,ii),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short \n daily mean T= 15 \n diurnal differnce= 9")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ddddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=31),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=31),fill="gray")+theme_bw()+geom_line(aes(xx,jj),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long \n daily mean T= 15 \n diurnal differnce= 2")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
noncovarying2<-ggarrange(aaaaa, bbbbb,ccccc, ddddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 

plot.list <- lapply(list(covarying,flat,noncovarying2,noncovarying), 
                    function(p) p + theme(plot.background = element_rect(color = "black")))
jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/designs.jpeg",width = 8,height=8,unit="in",res=300)
ggarrange(plotlist = plot.list,ncol=2,nrow=2,labels = c("a)","b)","c)","d)"))
dev.off()
grid.rect(width = 1, height = 0,.5, gp = gpar(lwd = 2, col = "black", fill = NA,hjust="left",vjust="topleft"))
grid.rect(width = 0,.5, height = 1, gp = gpar(lwd = 2, col = "black", fill = NA))
annotate_figure()

