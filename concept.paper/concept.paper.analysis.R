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
#mich.data$fruiting[mich.data$fruiting==19]<-7
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
#PhyloD <- phylo.d(d, binvar=pro) ###regular hysteranthy
#PhyloD
##functionalhysteranthy
PhyloPro2<-phylo.d(d,binvar=pro2)
PhyloPro2
#d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
#PhyloPro3<-phylo.d(d,binvar=pro3)
#PhyloPro3
###centering
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))

###prepare for modeling
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
mich.data$height10<-mich.data$heigh_height/10
###models:
mich5<-phyloglm(pro2~pol+height10+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5)

Mich5cent.funct<-phyloglm(pro2~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)

summary(Mich5cent.funct)

Mich5cent<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)
summary(Mich5cent)


Mich5cent.super<-phyloglm(pro3~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=50,full.matrix = TRUE)
#summary(Mich5cent.super)

###average predictive comparsons using uncentered data except for height because it wouldnt converge
beta<-coef(mich5)
hi<-1
lo<-0
###for [pollination syndrome]
beta[2]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$height10+beta[4]*mich.data$flo_time+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$height10+beta[4]*mich.data$flo_time+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)

print(mean(delta))

###for flowering
earl<-4
mid<-5
delta<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height10+beta[4]*earl+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height10+beta[4]*mid+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)

print(mean(delta))

#### restricted dataset (just early flowering)
median(mich.data$flo_time)


