rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)

GDH<-c(376,408,504,472,384,480,480,384)

photo<-c(11,15,15,11,11,11,15,15)
periodicities<-c("coupled","coupled","coupled","coupled","uncoupled","uncoupled","uncoupled","uncoupled") 

dat<-data.frame(GDH,photo,periodicities)
x<-c(10.8,10.8,15.2,15.2) 
y<-c(382,482,382,482)

xy<-data.frame(x,y)


ploty<-ggplot(dat,(aes(x=photo)))+geom_point(aes(y=GDH,color=periodicities),size=4)+xlim(8,18)+ylim(350,520)+theme_classic()+geom_polygon(aes(x=photo,y=GDH,color=periodicities,fill=periodicities),alpha=0.2)

plotty2<-ploty+geom_rect(dat,mapping=aes(xmin=10.9,xmax=15.1,ymin=382,ymax=482),fill=NA,color="black",size=0.2,linetype="dashed")

plotty2+geom_point(xy,mapping=aes(x,y),shape=8,size=4) 



12*24+12*16
12*20+12*12

24*15+16*9
24*11+16*13

20*15+12*9
20*11+12*13
))
