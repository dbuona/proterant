rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(gtable)
library(ggthemes)

GDH<-c(376,408,504,472,384,480,480,384)

photo<-c(11,15,15,11,11,11,15,15)
periodicities<-c("coupled","coupled","coupled","coupled","uncoupled","uncoupled","uncoupled","uncoupled") 

dat<-data.frame(GDH,photo,periodicities)
x<-c(10.8,10.8,15.2,15.2) 
y<-c(382,482,382,482)

xy<-data.frame(x,y)


ploty<-ggplot(dat,(aes(x=photo)))+geom_point(aes(y=GDH,color=periodicities),size=4)+xlim(8,18)+ylim(350,520)+theme_base()+geom_polygon(aes(x=photo,y=GDH,fill=periodicities),alpha=0.5)+scale_fill_manual(values=c("tomato1","royalblue2"))

plotty2<-ploty+geom_rect(dat,mapping=aes(xmin=10.9,xmax=15.1,ymin=382,ymax=482),fill=NA,color="black",size=0.2,linetype="dashed")

plotty2<-plotty2+geom_point(xy,mapping=aes(x,y),shape=8,size=4)+theme(legend.position = "left")


########
temp<-c(12,12,12,12)
light<-c(11,11,15,15)
treat<-c("warm/short","cool/short","warm/long","warm/long")
dawn<-data.frame(temp,light,treat)
dawn$dawn.offset<-dawn$temp-dawn$light
ggplot()


xx<-c(1:24)
yy<-rep(25,24)
zz<-rep(15,24)
aa<-c(20,20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20,20)
bb<-c(20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20)
cc<-c(10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10)
dd<-c(10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10)
ee<-c(20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20)
ff<-c(10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10)
gg<-c(20,20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20)
hh<-c(20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20,20) 
ii<-c(10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10)
jj<-c(10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10)

dat.simple<-data.frame(xx,yy,zz,aa,bb,cc,dd,ee,ff,gg,hh,ii)
a<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
b<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
c<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
d<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
flat<-ggarrange(a, b,c, d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

aaa<-ggplot(dat.simple,aes(xx,aa))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,aa),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,bb),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,cc),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,dd),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
covarying<-ggarrange(aaa, bbb,ccc, ddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 

aaaa<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
cccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
dddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
noncovarying<-ggarrange(aaaa, bbbb,cccc, dddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 


aaaaa<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,gg),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
bbbbb<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,hh),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ccccc<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=8,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=18,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,ii),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
ddddd<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=0,xmax=6,ymin=0,ymax=30),fill="gray")+geom_rect(dat.simple,mapping=aes(xmin=20,xmax=24,ymin=0,ymax=30),fill="gray")+theme_bw()+geom_line(aes(xx,jj),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))+theme(plot.title = element_text(face="italic", size=10))
noncovarying2<-ggarrange(aaaaa, bbbbb,ccccc, ddddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 
grid.arrange(flat,covarying,noncovarying,noncovarying2,ncol=2,nrow=2)
grid.rect(width = 1, height = 0,.5, gp = gpar(lwd = 2, col = "black", fill = NA,hjust="left",vjust="topleft"))
grid.rect(width = 0,.5, height = 1, gp = gpar(lwd = 2, col = "black", fill = NA))
annotate_figure()

