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
HF$name[HF$species=="ACRU"]<-"*A. rubrum"
HF$name[HF$species=="ACSA"]<-"*A. saccharrum"
HF$name[HF$species=="AMSP"]<-"Amelanchier spp."
HF$name[HF$species=="BEAL"]<-"*B. allegheniensis"
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
HF$name[HF$species=="QURU"]<-"*Q. rubra"
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
colnames(HF)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion (75%)","flower budburst","flower open","name","physiological offset","functional offset","intermediate offset")
HF2<-tidyr::gather(HF, phase,DOY,4:7)
HF2<-filter(HF2,!is.na(name))
jpeg("HarvardForest/HFmeans_expanded.jpeg",width = 8, height = 6, units = 'in', res=300)
ggplot(HF2,(aes(name,DOY)))+stat_summary(aes(shape=phase,color=phase))+scale_color_manual(values=c("darkgray","darkgray","black","black"))+scale_shape_manual(values=c(2,1,17,16))+theme_linedraw()+ylab("Day of Year")+xlab(NULL)+theme(axis.text.x = element_text(angle = 300,hjust=0,face="italic"))
dev.off()


########Analysis 3 Quercus rubra at HF.........
QURU<-filter(HF,species=="QURU")
QURU<-filter(QURU,tree.id!="QURU-02")
colnames(QURU)
QURU$FLS<-ifelse(QURU$"physiological offset"<0,"seranthous","hysteranthous")
#colnames(QURU)<-c("year","tree.id","species","leaf_budburst","leaf_expansion(75%)","flower_budburst","flower_open","physiological_offset","functional_offset","name","FLS")
QURU<-drop_na(QURU)
QURU$tree.id[QURU$tree.id=="QURU-01"]<-"1"
QURU$tree.id[QURU$tree.id=="QURU-03"]<-"3"
QURU$tree.id[QURU$tree.id=="QURU-04"]<-"4"

colnames(QURU)[6]<-"flower.budburst"
colnames(QURU)[4]<-"leaf.budburst"
jpeg("HarvardForest/HF_Q_ru_interannual.jpeg",width = 8, height = 6, units = 'in', res=350)
pd2<-position_dodge2(0.7)

ggplot(QURU,aes(year,flower.budburst))+geom_point(aes(year,flower.budburst,color=tree.id,group = row.names(QURU),shape="flower.budburst") ,position=pd2)+
  geom_point(aes(year,leaf.budburst,color=tree.id,group = row.names(QURU),shape="leaf.budburst"),position=pd2)+
geom_linerange(aes(x=year,ymin=flower.budburst,ymax=leaf.budburst, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd2)+
  theme_linedraw()+labs(y = "Day of year",color= "Tree I.D.")+
scale_color_brewer(palette="Set1")+
  scale_shape_manual(name="Phenophase",values=c(2,17), label=c("flower budburst","leaf budburst"))
dev.off()

