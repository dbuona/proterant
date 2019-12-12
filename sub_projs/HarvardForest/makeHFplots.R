rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_projs")


library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")
library(tidyr)


##read in the data
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)

###give them species names
HF$name[HF$species=="ACPE"]<-"A. pensylvanicum"
HF$name[HF$species=="ACRU"]<-"A. rubrum"
HF$name[HF$species=="ACSA"]<-"*A. saccharrum"
HF$name[HF$species=="AMSP"]<-"Amelanchier spp."
HF$name[HF$species=="BEAL"]<-"B. allegheniensis"
HF$name[HF$species=="BELE"]<-"*B. lenta"
HF$name[HF$species=="BEPA"]<-"*B. papyrifera"
HF$name[HF$species=="BEPO"]<-"*B. populifolia"
HF$name[HF$species=="FAGR"]<-"*F. grandifolia"
HF$name[HF$species=="FRAM"]<-"*F. americana"
HF$name[HF$species=="ILVE"]<-"I. verticillata"
HF$name[HF$species=="KAAN"]<-"K. angustifolia"
HF$name[HF$species=="KALA"]<-"K. latifolia"
HF$name[HF$species=="NEMU"]<-"I. mucronata"
HF$name[HF$species=="NYSY"]<-"N. sylvatica"
HF$name[HF$species=="POTR"]<-"*P. tremuloides"
HF$name[HF$species=="PRSE"]<-"P. serotina"
HF$name[HF$species=="QURU"]<-"Q. rubra"
HF$name[HF$species=="QUVE"]<-"*Q. velutina"
HF$name[HF$species=="VACO"]<-"V. corymbosum"
HF$name[HF$species=="SAPU"]<-"S. racemosa"
HF$name[HF$species=="VICA"]<-"V. cassinoides"
HF$name[HF$species=="VIAL"]<-"V. latinoides"
HF$name[HF$species=="COAL"]<-"C. alternifolia"
HF$phys.fls<-HF$bb.jd-HF$fbb.jd
HF$funct.fls<-HF$l75.jd-HF$fopn.jd
HF$inter.fls<-HF$bb.jd-HF$fopn.jd

colnames(HF)
colnames(HF)<-c("year" , "tree.id", "species","Lbb","L75","Fbb","Fopn","name","physiological offset","functional offset","intermediate offset")
HF2<-tidyr::gather(HF, Phenophase,DOY,4:7)
HF2<-filter(HF2,!is.na(name))
HFshort<-dplyr::filter(HF2,name %in% c("A. rubrum","Amelanchier spp."))

#jpeg("HarvardForest/HFmeans_expanded.jpeg",width = 8, height = 6, units = 'in', res=300)
a<-ggplot(HFshort,(aes(name,DOY)))+stat_summary(aes(shape=Phenophase,color=Phenophase),size=1)+scale_color_manual(values=c("darkgray","darkgray","black","black"))+scale_shape_manual(values=c(2,1,17,16))+theme_linedraw()+ylab("Day of Year")+xlab(NULL)+theme(axis.text.x = element_text(face="italic",angle=300,hjust = 0.1))+
  labs(title = NULL,tag="a)")+ theme(legend.position = "bottom",legend.box = "vertical")


####does the graph work for these species
shorty<-filter(HF,name %in% c("Q. rubra"))
#shorty<-filter(shorty,year<2002)
colnames(shorty)
shorty$FLS<-ifelse(shorty$"physiological offset"<0,"seranthous","hysteranthous")
#colnames(shorty)<-c("year","tree.id","species","leaf_budburst","leaf_expansion(75%)","flower_budburst","flower_open","physiological_offset","functional_offset","name","FLS")
shorty<-drop_na(shorty)
shorty$tree.id[shorty$tree.id=="QURU-01"]<-"1"
shorty$tree.id[shorty$tree.id=="QURU-02"]<-"2"
shorty$tree.id[shorty$tree.id=="QURU-03"]<-"3"
shorty$tree.id[shorty$tree.id=="QURU-04"]<-"4"

colnames(shorty)[6]<-"flower.budburst"
colnames(shorty)[4]<-"leaf.budburst"
#jpeg("HarvardForest/HF_Q_ru_interannual.jpeg",width = 8, height = 6, units = 'in', res=350)
pd2<-position_dodge2(0.6)

b<-ggplot(shorty,aes(year,flower.budburst))+geom_point(aes(year,flower.budburst,color=tree.id,group = row.names(shorty),shape="flower.budburst") ,position=pd2)+
  geom_point(aes(year,leaf.budburst,color=tree.id,group = row.names(shorty),shape="leaf.budburst"),position=pd2)+
  geom_linerange(aes(x=year,ymin=flower.budburst,ymax=leaf.budburst, linetype=FLS,color=tree.id,group = row.names(shorty)),position=pd2)+
  theme_linedraw()+labs(y = "Day of year",color= "Tree I.D.")+
  scale_color_brewer(palette="Set1")+
  scale_shape_manual(name="Phenophase",values=c(2,17), label=c("flower budburst","leaf budburst"))+
  labs(title = NULL,tag="b)")


jpeg("HarvardForest/FLS_viz.jpeg",width = 10, height = 4.4, units = 'in', res=200)
ggpubr::ggarrange(a,b,widths=c(.7,1))
dev.off()


#jpeg("HarvardForest/HF_Q_ru_interannual.jpeg",width = 8, height = 6, units = 'in', res=350)
