###new hysteranthy figure
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

setwd("~/Documents/git/proterant/sub_projs/")
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(tibble)
library(raster)
library(lme4)
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
library("ggrepel")


frax50<-read.csv("PEP725/input/frax50year.csv")




alnumans<- frax50 %>% group_by(s_id,lat,lon) %>% summarise(meanFLS=mean(FLS),sdFLS=sd(FLS)) 
alnumans$meanFLS<-(floor(alnumans$meanFLS))
alnumans$sdFLS<-(floor(alnumans$sdFLS))
alnumans$hyst<-ifelse(alnumans$meanFLS>0,"hyst","ser")
#alnumans$sdFLS <- paste0("(",alnumans$sdFLS , ")")
alnumans$FLS <- paste0(alnumans$meanFLS,alnumans$sdFLS)
randomRows <- function(df,n){
  return(df[sample(nrow(df),n),])
}
shortal<-randomRows(alnumans,20)

world <- ne_countries(scale = "medium", returnclass = "sf")




a<-ggplot(data = world) + geom_sf(fill="white",color="gray")+geom_text_repel(data = shortal, aes(x = lon, y = lat, label = FLS),size=3, 
                                                 nudge_x = c(.5, -.5, 1, 2, -1), nudge_y = c(0.25, -0.25, 0.5, 0.5, -0.5)) +
  geom_point(data = alnumans, aes(x = lon, y = lat), size = .5)+
  coord_sf(xlim = c(5.5, 14.5), ylim = c(47, 55), expand = TRUE)+ggthemes::theme_base()

setEPS()
postscript("fraxmeans.eps",width = 6, height = 4)
ggplot(alnumans,aes(meanFLS))+geom_histogram(binwidth = 1,color="black",fill="lightgray")+ggthemes::theme_base(base_size = 11)
dev.off()
setEPS()
postscript("fraxSD.eps",width = 6, height = 4)
ggplot(alnumans,aes(sdFLS))+geom_histogram(binwidth = 1,color="black",fill="lightgray")+ggthemes::theme_base(base_size = 11)
dev.off()



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
HF$name[HF$species=="BEPA"]<-"B. papyrifera"
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
HF$name[HF$species=="QUVE"]<-"Q. velutina"
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


hyst<-c("A. rubrum", "F. americana")
syn<-c("A. pensylvanicum","Q. rubra")
ser<-c("N. sylvatica","P. serotina")

HFshort$FLS<-NA
HFshort$FLS[which(HFshort$name %in% hyst)] <- "hysteranthous"
HFshort$FLS[which(HFshort$name %in% syn)] <- "synanthous"
HFshort$FLS[which(HFshort$name %in% ser)] <- "seranthous"

HFshort$FLS = factor(HFshort$FLS, levels=c("hysteranthous","synanthous","seranthous"))

b<-HFshort%>%
  mutate(name = factor(name, levels=c("A. rubrum","F. americana", "A. pensylvanicum","Q. rubra","N. sylvatica","P. serotina"))) %>%
  ggplot((aes(name,DOY)))+stat_summary(aes(shape=Phenophase,color=Phenophase))+scale_color_manual(values=c("darkgray","darkgray","black","black"))+scale_shape_manual(values=c(2,1,17,16))+ggthemes::theme_base(base_size = 11)+ylab("Day of Year")+xlab(NULL)+theme(axis.text.x = element_text(face="italic",angle=340,vjust = -0.4,hjust=0.3))+
  labs(tag="a)")+facet_wrap(~FLS,scale="free_x",strip.position  ="top")



Arub<-filter(HFshort,name=="A. rubrum")
colnames(Arub)[5]<-"phys"
colnames(Arub)[6]<-"funct"
colnames(Arub)[7]<-"inter"
c<-ggplot(data=Arub,aes(x=tree.id,phys))+geom_boxplot()+ggthemes::theme_base(base_size = 10)+labs(tag="c)")+scale_x_discrete(name="Indiviudals",label=c(1,2,3,4,5))
?scale_x_discrete()

shorty<-filter(HF,tree.id %in% c("QURU-04","QURU-03","QURU-01"))
table(shorty$tree.id)
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
pd2<-position_dodge2(0.8)
colourCount = length(unique(shorty$tree.id))





d<-ggplot(shorty,aes(year,flower.budburst))+geom_point(aes(year,flower.budburst,color=tree.id,group = row.names(shorty),shape="flower.budburst") ,position=pd2,size=2)+
  geom_point(aes(year,leaf.budburst,color=tree.id,group = row.names(shorty),shape="leaf.budburst"),position=pd2,size=2)+
  geom_linerange(aes(x=year,ymin=flower.budburst,ymax=leaf.budburst,color=tree.id, linetype=FLS,group = row.names(shorty)),position=pd2)+
  ggthemes::theme_base(base_size = 10)+labs(y = "Day of year",color= "Tree I.D.")+guides(color = FALSE)+scale_color_brewer(palette = "Set1",type="quan")+
  scale_shape_manual(name="Phenophase",values=c(2,17), label=c("flower budburst","leaf budburst"))+
  labs(tag="b)")+theme(strip.text.x = element_text(face="italic"))#+
  theme(axis.text.x = element_text(angle=300,hjust = 0.1))

  setEPS()
  postscript("intraspecificplots.eps",width = 8, height = 10)
ggpubr::ggarrange(b,d,nrow=2,ncol=1,widths=c(1,1))
dev.off()
ggpubr::ggarrange(goo,d,nrow=2,heights = c(1,1))

setEPS()
postscript("popmap.eps",width = 8, height = 10)
a
dev.off()
