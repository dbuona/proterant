rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_projs")


library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")
library(tidyr)
library(RColorBrewer)
library(grid)

##read in the data
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)

###give them species names
HF$name[HF$species=="ACPE"]<-"A. pensylvanicum"
HF$name[HF$species=="ACRU"]<-"A. rubrum"
HF$name[HF$species=="ACSA"]<-"A. saccharrum"
HF$name[HF$species=="AMSP"]<-"Amelanchier spp."
HF$name[HF$species=="BEAL"]<-"B. allegheniensis"
HF$name[HF$species=="BELE"]<-"B. lenta"
HF$name[HF$species=="BEPA"]<-"*B. papyrifera"
HF$name[HF$species=="BEPO"]<-"B. populifolia"
HF$name[HF$species=="FAGR"]<-"F. grandifolia"
HF$name[HF$species=="FRAM"]<-"F. americana"
HF$name[HF$species=="ILVE"]<-"I. verticillata"
HF$name[HF$species=="KAAN"]<-"K. angustifolia"
HF$name[HF$species=="KALA"]<-"K. latifolia"
HF$name[HF$species=="NEMU"]<-"I. mucronata"
HF$name[HF$species=="NYSY"]<-"N. sylvatica"
HF$name[HF$species=="POTR"]<-"P. tremuloides"
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
colnames(HF)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion","flower budburst","flowers open","name","physiological offset","functional offset","intermediate offset")
HF2<-tidyr::gather(HF, Phenophase,DOY,4:7)
HF2<-filter(HF2,!is.na(name))
HFshort<-dplyr::filter(HF2,name %in% c("A. rubrum","F. americana","Q. rubra","A. pensylvanicum","P. serotina","N. sylvatica"))


#jpeg("HarvardForest/HFmeans_expanded.jpeg",width = 8, height = 6, units = 'in', res=300)
dater<-data_frame(name=rep(unique(HFshort$name),each=2),FLS=rep(c("first","second"),6),Phenophase=rep(c("flower","leaf"),6))

#a<-ggplot()+ annotate("text", x = 5, y = 25, label = "Flowers before leaves",size=5,fontface=2)+
 # annotate("text", x = 5, y = 24.9, label = "A. rubrum",size=3,fontface=3 )+
#  annotate("text", x = 5, y = 24.85, label = "F. americana",size=3,fontface=3 )+
##  annotate("text", x = 5, y = 24.5, label = "Flowers with leaves",size=5,fontface=2)+
#  annotate("text", x = 5, y = 24.4, label = "A. pensylvanicum",size=3,fontface=3 )+
#  annotate("text", x = 5, y = 24.35, label = "Q. rubra",size=3,fontface=3 )+
#  annotate("text", x = 5, y = 24, label = "Flowers after leaves",size=5,fontface=2)+
#annotate("text", x = 5, y = 23.9, label = "N. sylvatica" ,size=3,fontface=3 )+
#  annotate("text", x = 5, y = 23.85, label = "P. serotina",size=3,fontface=3 )+
#  theme_bw()+theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())+
#labs(title = "Current definitions",tag="a)")

hyst<-c("A. rubrum", "F. americana")
syn<-c("A. pensylvanicum","Q. rubra")
ser<-c("N. sylvatica","P. serotina")

HFshort$FLS<-NA
HFshort$FLS[which(HFshort$name %in% hyst)] <- "hysteranthous"
HFshort$FLS[which(HFshort$name %in% syn)] <- "synanthous"
HFshort$FLS[which(HFshort$name %in% ser)] <- "seranthous"

HFshort$FLS = factor(HFshort$FLS, levels=c("hysteranthous","synanthous","seranthous"))

a<-HFshort%>%
  mutate(name = factor(name, levels=c("A. rubrum","F. americana", "A. pensylvanicum","Q. rubra","N. sylvatica","P. serotina"))) %>%
ggplot((aes(name,DOY)))+stat_summary(aes(shape=Phenophase,color=Phenophase))+scale_color_manual(values=c("darkgray","darkgray","black","black"))+scale_shape_manual(values=c(2,1,17,16))+theme_linedraw()+ylab("Day of Year")+xlab(NULL)+theme(axis.text.x = element_text(face="italic",angle=300,hjust = 0.1))+
  labs(title = "Quantitative phenology in the field",tag="a)")+facet_wrap(~FLS,scale="free_x",strip.position  ="top")



# theme(legend.position = "bottom",legend.box = "vertical")+

#AA<-ggpubr::ggarrange(a,b,widths=c(1.1,2))
####does the graph work for these species
shorty<-filter(HF,name %in% c("F. americana","Q. rubra"))
#shorty<-filter(shorty,year<2002)
colnames(shorty)
shorty$FLS<-ifelse(shorty$"physiological offset"<0,"seranthous","hysteranthous")
#colnames(shorty)<-c("year","tree.id","species","leaf_budburst","leaf_expansion(75%)","flower_budburst","flower_open","physiological_offset","functional_offset","name","FLS")
shorty<-drop_na(shorty)
#shorty$tree.id[shorty$tree.id=="QURU-01"]<-"1"
#shorty$tree.id[shorty$tree.id=="QURU-02"]<-"2"
#shorty$tree.id[shorty$tree.id=="QURU-03"]<-"3"
#shorty$tree.id[shorty$tree.id=="QURU-04"]<-"4"

colnames(shorty)[6]<-"flower.budburst"
colnames(shorty)[4]<-"leaf.budburst"
#jpeg("HarvardForest/HF_Q_ru_interannual.jpeg",width = 8, height = 6, units = 'in', res=350)
pd2<-position_dodge2(0.6)
colourCount = length(unique(shorty$tree.id))


b<-ggplot(shorty,aes(year,flower.budburst))+geom_point(aes(year,flower.budburst,color=tree.id,group = row.names(shorty),shape="flower.budburst") ,position=pd2)+
  geom_point(aes(year,leaf.budburst,color=tree.id,group = row.names(shorty),shape="leaf.budburst"),position=pd2)+
  geom_linerange(aes(x=year,ymin=flower.budburst,ymax=leaf.budburst, linetype=FLS,color=tree.id,group = row.names(shorty)),position=pd2)+
  theme_linedraw()+labs(y = "Day of year",color= "Tree I.D.")+guides(color = FALSE)+
  scale_shape_manual(name="Phenophase",values=c(2,17), label=c("flower budburst","leaf budburst"))+
  labs(title = "Intra-specific variation",tag="b)")+facet_wrap(~name)+theme(strip.text.x = element_text(face="italic"))+
  theme(axis.text.x = element_text(angle=300,hjust = 0.1))



ggpubr::ggarrange(a,b,ncol=2,widths=c(1.5,2))

jpeg("HarvardForest/FLS_viz.jpeg",width = 10, height = 8, units = 'in',res=350)
ggpubr::ggarrange(a,b,nrow=2,heights=c(2.2,2))

dev.off()


#jpeg("HarvardForest/HF_Q_ru_interannual.jpeg",width = 8, height = 6, units = 'in', res=350)
