####concept paper code: THis is based off hyst_final_analysis.R. but only what we need to the manuscript
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")
###cleaning:
#clean av fruit time
mich.data$av_fruit_time[mich.data$av_fruit_time=="persistant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="persitant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="unreported"]<-9    
mich.data$av_fruit_time<-as.numeric(mich.data$av_fruit_time)

mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9                                      
mich.data$fruiting<-as.numeric(mich.data$fruiting)

mich.data["pro3"]<-NA
mich.data$pro3[mich.data$Phen.sequence == "pro"] <- 1
mich.data$pro3[mich.data$Phen.sequence == "pro/syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn/ser"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "ser"] <- 0 
mich.data$pro3[mich.data$Phen.sequence== "hyst"] <- 0


###phylo.D
set.seed(122)
mich.tree$node.label<-NULL
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloD <- phylo.d(d, binvar=pro) ###regular hysteranthy
PhyloD
##functionalhysteranthy
PhyloPro2<-phylo.d(d,binvar=pro2)
PhyloPro2
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloPro3<-phylo.d(d,binvar=pro3)
PhyloPro3
###centering
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))

###prepare for modeling
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")

###models:
mich5<-phyloglm(pro~pol+heigh_height+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)

mich5cent<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)

Mich5cent.funct<-phyloglm(pro2~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)

Mich5cent.super<-phyloglm(pro3~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)
###effect plots
bootest<-as.data.frame(Mich5cent.funct$coefficients)
bootconf<-as.data.frame(Mich5cent.funct$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="alpha")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
bootmich$trait[bootmich$trait=="shade_bin"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="fruit_cent"]<-"fruit timing"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

functplot<-ggplot(bootmich,aes(estimate,trait))+geom_point()+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+ggtitle("Functional Hysteranthy")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+theme(plot.title = element_text(hjust = 0.5))+guides(color="none")
####now phys
bootest<-as.data.frame(Mich5cent.super$coefficients)
bootconf<-as.data.frame(Mich5cent.super$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="alpha")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
bootmich$trait[bootmich$trait=="shade_bin"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="fruit_cent"]<-"fruit timing"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

physplot<-ggplot(bootmich,aes(estimate,trait))+geom_point()+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+ggtitle("Physiological Hysteranthy")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+theme(plot.title = element_text(hjust = 0.5))+guides(color="none")

