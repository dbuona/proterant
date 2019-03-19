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


ploty<-ggplot(dat,(aes(x=photo)))+geom_point(aes(y=GDH,color=periodicities),size=4)+xlim(8,18)+ylim(350,520)+theme_base()+geom_polygon(aes(x=photo,y=GDH,color=periodicities,fill=periodicities),alpha=0.5) +scale_color_manual(values=c("rosybrown","gray10"))+scale_fill_manual(values=c("rosybrown","gray10"))

plotty2<-ploty+geom_rect(dat,mapping=aes(xmin=10.9,xmax=15.1,ymin=382,ymax=482),fill=NA,color="black",size=0.2,linetype="dashed")

plotty2<-plotty2+geom_point(xy,mapping=aes(x,y),shape=8,size=4)+theme(legend.position = "left")


xx<-c(1:24)
yy<-rep(25,24)
zz<-rep(15,24)
aa<-c(30,30,30,30,30,30,30,30,30,30,20,20,20,20,20,20,20,20,20,20,20,20,20,20)
bb<-c(30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20,20,20,20)
cc<-c(20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
dd<-c(20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,10)
ee<-c(30,30,30,30,30,30,30,30,30,30,30,30,20,20,20,20,20,20,20,20,20,20,20,20)
ff<-c(20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,10,10,10,10,10)

dat.simple<-data.frame(xx,yy,zz,aa,bb,cc,dd,ee,ff)
a<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=10,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
b<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,yy),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
c<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=10,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
d<-ggplot(dat.simple,aes(xx,yy))+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,zz),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
flat<-ggarrange(a, b,c, d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
?ggarrange()

aaa<-ggplot(dat.simple,aes(xx,aa))+geom_rect(dat.simple,mapping=aes(xmin=10,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,aa),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
bbb<-ggplot(dat.simple,aes(xx,bb))+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,bb),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
ccc<-ggplot(dat.simple,aes(xx,cc))+geom_rect(dat.simple,mapping=aes(xmin=10,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,cc),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
ddd<-ggplot(dat.simple,aes(xx,dd))+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,dd),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
covarying<-ggarrange(aaa, bbb,ccc, ddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 

aaaa<-ggplot(dat.simple,aes(xx,ee))+geom_rect(dat.simple,mapping=aes(xmin=10,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
bbbb<-ggplot(dat.simple,aes(xx,ee))+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,ee),color="red")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("warm/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
cccc<-ggplot(dat.simple,aes(xx,ff))+geom_rect(dat.simple,mapping=aes(xmin=10,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/short")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
dddd<-ggplot(dat.simple,aes(xx,ff))+geom_rect(dat.simple,mapping=aes(xmin=16,xmax=24,ymin=0,ymax=30),fill="gray")+theme_base()+geom_line(aes(xx,ff),color="blue")+scale_x_continuous(breaks=c(4,8,12,16,20,24))+scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+ labs(y = "temperature",x = "hours")+ggtitle("cool/long")+theme(plot.title=element_text( hjust=0.5, vjust=0.5))
noncovarying<-ggarrange(aaaa, bbbb,cccc, dddd, ncol=2, nrow=2, common.legend = TRUE, legend="bottom") 

grid.arrange(flat,covarying,noncovarying,ncol=3)
grid.rect(width = .5, height = .5, gp = gpar(lwd = 2, col = "black", fill = NA,hjust="left",vjust="top"))
grid.rect(width = 1, height = 1, gp = gpar(lwd = 2, col = "black", fill = NA))


