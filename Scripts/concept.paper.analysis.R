####c### This is the final analysis file for hysteranthy anaylsis on MTSV as of 3/28/18.
###major update on 11.27.18
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
library("randomForest")
library(car)
library(ggstance)
library(broom)
library(brms)
library(rstan)

#if you dont want to run the model: 
load("RData/zarchival/hystmodels.RData")
#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
drought.dat<-read.csv("..//Data/USDA_traitfor_MTSV.csv",header=TRUE)
source("..//Scripts/extract_coefs.R")
###make a column for seed development time
mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time
###one more cleaninging tax
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)
###make the tree work 
mich.tree$node.label<-NULL

###RESPONSE VARIABLE KEY
#pro2<-hysteranthy= before, before/with and with
#pro3<- hysteranthy=before only
###phylo.D###### calculate phylo d.

set.seed(122)
######phylo signals#########################################################
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
Wind<-phylo.d(d, binvar=pol) #-0.48 more phylogenetically conserved that 
Wind
##functionalhysteranthy
PhyloPro2<-phylo.d(d,binvar=pro2)
PhyloPro2
plot(PhyloPro2)

PhyloPro3<-phylo.d(d,binvar=pro3)
PhyloPro3
plot(PhyloPro3)
###phlosignal for continuous trait
phylosig(mich.tree, mich.data$flo_time, method="lambda", test=TRUE, nsim=999,se=NULL)
###############################################################################################
######### later we'll be running model with drought tolerance, so add drought tolerance data
mich.data<-left_join(mich.data,drought.dat,by="name")

######prune the tree for drought modeling
mich.data<-filter(mich.data,!is.na(min._precip))

####prune tree to match reduced dataset
names.intree<-mich.tree$tip.label
namelist<-unique(mich.data$name)
to.prune<-which(!names.intree%in%namelist)
mich.tree.droughtprune<-drop.tip(mich.tree,to.prune)
mytree.names<-mich.tree.droughtprune$tip.label
### the smaller new tree is called mich.tree.droughtprune

####recenter here after you prune the list
###Rescale predictors this makes it so you can compare binary to continous data
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))
mich.data$dev_time_cent<-(mich.data$dev.time-mean(mich.data$dev.time))/(2*sd(mich.data$dev.time))
mich.data$tol_cent<-(mich.data$shade_bin-mean(mich.data$shade_bin))/(2*sd(mich.data$shade_bin))
mich.data$precip_cent<-(mich.data$min._precip-mean(mich.data$min._precip))/(2*sd(mich.data$min._precip))


##############################################
###Now do all you did with silvics
#######################################################
silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")
silv.USDA<-read.csv("silv.USDA.csv")
silv.tree$node.label<-NULL
silv.data<-left_join(silv.data,silv.USDA) ###maybe you should make a drought species list specifically for silvics to lose less species

####Silvics cleaning########
###fruiting
silv.data$fruiting<-NA
silv.data$fruiting<-silv.data$av_fruit_time
#silv.data$fruiting[silv.data$fruiting==21]<-9

###functional hysteranthy
silv.data["pro2"]<-NA
silv.data$pro2[silv.data$silvic_phen_seq== "pro"] <- 1
silv.data$pro2[silv.data$silvic_phen_seq== "pro/syn"] <- 1
silv.data$pro2[silv.data$silvic_phen_seq== "syn"] <- 1
silv.data$pro2[silv.data$silvic_phen_seq== "syn/ser"] <- 0
silv.data$pro2[silv.data$silvic_phen_seq== "ser"] <- 0 
silv.data$pro2[silv.data$silvic_phen_seq== "hyst"] <- 0
silv.data$pro2[silv.data$name == "Quercus_laurifolia"] <- 1

###super bio hysteranthy for silvics
silv.data["pro3"]<-NA
silv.data$pro3[silv.data$silvic_phen_seq== "pro"] <- 1
silv.data$pro3[silv.data$silvic_phen_seq== "pro/syn"] <- 0
silv.data$pro3[silv.data$silvic_phen_seq== "syn"] <- 0
silv.data$pro3[silv.data$silvic_phen_seq== "syn/ser"] <- 0
silv.data$pro3[silv.data$silvic_phen_seq== "ser"] <- 0 
silv.data$pro3[silv.data$silvic_phen_seq== "hyst"] <- 0
silv.data$pro3[silv.data$name == "Quercus_laurifolia"] <- 1

silv.data<-filter(silv.data,!is.na(min._precip)) 
###silv rescaling
silv.data$height_cent<-(silv.data$height-mean(silv.data$height))/(2*sd(silv.data$height))
silv.data$fruit_cent<-(silv.data$fruiting-mean(silv.data$fruiting))/(2*sd(silv.data$fruiting))
silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))
silv.data$pol_cent<-(silv.data$pol-mean(silv.data$pol))/(2*sd(silv.data$pol))
silv.data$tol_cent<-(silv.data$shade_bin-mean(silv.data$shade_bin))/(2*sd(silv.data$shade_bin))
silv.data$dev.time<-silv.data$fruiting-silv.data$flower_time
silv.data$dev_time_cent<-(silv.data$dev.time-mean(silv.data$dev.time))/(2*sd(silv.data$dev.time))
silv.data$precip_cent<-(silv.data$min._precip-mean(silv.data$min._precip))/(2*sd(silv.data$min._precip)) 

####prune silvics tree to match reduced dataset
names.intree<-silv.tree$tip.label
namelist<-unique(silv.data$name)
to.prune<-which(!names.intree%in%namelist)
silv.tree.droughtprune<-drop.tip(silv.tree,to.prune)
mytree.names<-silv.tree.droughtprune$tip.label
setdiff(namelist,mytree.names) 
intersect(namelist,mytree.names)

#####phylogenetic signals####################################### for silvics
e<-comparative.data(silv.tree,silv.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloSilv<-phylo.d(e,binvar=pro)
PhyloSilv
PhyloSilv2<-phylo.d(e,binvar=pro2)
PhyloSilv2
plot(PhyloSilv2)
PhyloSilv3<-phylo.d(e,binvar=pro3)
PhyloSilv3
plot(PhyloSilv3)
#############################################################################################
#### interaction models 
##########################################################################################
######### phyloglm requires species names to be in rownames
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<- silv.data %>% remove_rownames %>% column_to_rownames(var="name")
###models: #number of bootstraps from Wilcox, R. R. (2010). Fundamentals of modern statistical methods: Substantially improving power and accuracy. Springer.

z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)


z.phys.drought<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

z.funct.drought.silvics<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                  start.beta=NULL, start.alpha=NULL,
                                  boot=599,full.matrix = TRUE)

z.phys.drought.silvics<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                 start.beta=NULL, start.alpha=NULL,
                                 boot=599,full.matrix = TRUE)
#======clean for plotting
mich.funct.wint.dat<-full_join(extract_coefs(z.funct.drought),extract_CIs(z.funct.drought),by="trait")
colnames(mich.funct.wint.dat)<-c("trait","estimate","low","high")
mich.funct.wint.dat$class<-"functional-MTSV"

mich.phys.wint.dat<-full_join(extract_coefs(z.phys.drought),extract_CIs(z.phys.drought),by="trait")
colnames(mich.phys.wint.dat)<-c("trait","estimate","low","high")
mich.phys.wint.dat$class<-"physiological-MTSV"

silv.funct.wint.dat<-full_join(extract_coefs(z.funct.drought.silvics),extract_CIs(z.funct.drought.silvics),by="trait")
colnames(silv.funct.wint.dat)<-c("trait","estimate","low","high")
silv.funct.wint.dat$class<-"functional-USFS"

silv.phys.wint.dat<-full_join(extract_coefs(z.phys.drought.silvics),extract_CIs(z.phys.drought.silvics),by="trait")
colnames(silv.phys.wint.dat)<-c("trait","estimate","low","high")
silv.phys.wint.dat$class<-"physiological-USFS"

michigan.wint<-rbind(mich.phys.wint.dat,mich.funct.wint.dat)
USFS.wint<-rbind(silv.phys.wint.dat,silv.funct.wint.dat)

#============== cleaning
comps<-rbind(michigan.wint,USFS.wint)
comps$category<-NA
comps$category[which(comps$class=="physiological-USFS")] <- "physiological"
comps$category[which(comps$class=="physiological-MTSV")] <- "physiological"
comps$category[which(comps$class=="functional-USFS")] <- "functional"
comps$category[which(comps$class=="functional-MTSV")] <- "functional"

comps$data<-NA
comps$data[which(comps$class=="physiological-USFS")] <- "USFS"
comps$data[which(comps$class=="physiological-MTSV")] <- "MTSV"
comps$data[which(comps$class=="functional-USFS")] <- "USFS"
comps$data[which(comps$class=="functional-MTSV")] <- "MTSV"

###change the variable names
comps$trait[which(comps$trait=="pol_cent")] <- "main effect: pollination syndrome"
comps$trait[which(comps$trait=="flo_cent")] <- "main effect: flowering time"
comps$trait[which(comps$trait=="precip_cent")] <- "main effect: minimum precipitation"
comps$trait[which(comps$trait=="flo_cent:precip_cent")] <- "interaction: flowering x precip."
comps$trait[which(comps$trait=="pol_cent:precip_cent")] <- "interaction: pollination x precip."
comps$trait[which(comps$trait=="pol_cent:flo_cent")] <- "interaction: pollination x flowering"

pd=position_dodgev(height=0.4)
plotty3<-ggplot(comps,aes(estimate,trait))+geom_point(size=2.5,aes(color=category,shape=data),position=pd)+geom_errorbarh(position=pd,width=0.4,aes(xmin=low,xmax=high,color=category,shape=data))+geom_vline(aes(xintercept=0))+theme_bw()+scale_color_manual(values=c("orchid4", "springgreen4"))
plotty3
###does the plot look better with the two separet
comps.MTSV<-filter(comps,data=="MTSV")
comps.USFS<-filter(comps,data=="USFS")
pd=position_dodgev(height=0.4)
plotty3a<-ggplot(comps.MTSV,aes(estimate,trait))+geom_point(size=2.5,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_bw()+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)
plotty3b<-ggplot(comps.USFS,aes(estimate,trait))+geom_point(size=2.5,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_bw()+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)
grid.arrange(plotty3a,plotty3b,nrow=1)

#models=====================no interaction=========================================================================================
z.funct.drought.noint<-phyloglm(pro2~pol_cent+flo_cent+precip_cent,mich.data, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

summary(z.funct.drought.noint)


####March 13, 2019 Start testing stan models
z.no.phylo.freq<-glm(pro2~pol_cent+flo_cent+precip_cent,mich.data,family=binomial(link="logit"))
summary(z.no.phylo.freq)
colnames(mich.data)
no.phylo.freq<-glm(pro2~pol+flo_time+min._precip,mich.data,family=binomial(link="logit"))
summary(no.phylo.freq)

?glm()
####Bayesian
datalist<- with(mich.data, 
                list(y=pro2,
                     pol=pol_cent,
                     flotime=flo_cent,
                     minP=precip_cent,
                     N = nrow(mich.data)
                ))

datalist2<- with(mich.data, 
                list(y=pro2,
                     pol=pol,
                     flotime=flo_time,
                     minP=min._precip,
                     N = nrow(mich.data)
                ))

z.no.phylo.bayes<- stan('..//stan/binary_stan_nophylo.stan', data = datalist,
             iter = 5000, warmup=3500) 

sum.bayes<-summary(z.no.phylo.bayes)$summary
sum.bayes[c("alpha","b_pol","b_flotime","b_minP"),]


no.phylo.bayes<- stan('..//stan/binary_stan_nophylo.stan', data = datalist2,
                        iter = 5000, warmup=3500)

sum.bayes2<-summary(no.phylo.bayes)$summary
sum.bayes2[c("alpha","b_pol","b_flotime","b_minP"),]

z.phys.drought.noint<-phyloglm(pro3~pol_cent+flo_cent+precip_cent,mich.data, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

z.funct.drought.silvics.noint<-phyloglm(pro2~pol_cent+flo_cent+precip_cent,silv.data, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                  start.beta=NULL, start.alpha=NULL,
                                  boot=599,full.matrix = TRUE)

z.phys.drought.silvics.noint<-phyloglm(pro3~pol_cent+flo_cent+precip_cent,silv.data, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                 start.beta=NULL, start.alpha=NULL,
                                 boot=599,full.matrix = TRUE)

#=============cleaning====================================================
mich.funct.noint.dat<-full_join(extract_coefs(z.funct.drought.noint),extract_CIs(z.funct.drought.noint),by="trait")
colnames(mich.funct.noint.dat)<-c("trait","estimate","low","high")
mich.funct.noint.dat$class<-"functional-MTSV"

mich.phys.noint.dat<-full_join(extract_coefs(z.phys.drought.noint),extract_CIs(z.phys.drought.noint),by="trait")
colnames(mich.phys.noint.dat)<-c("trait","estimate","low","high")
mich.phys.noint.dat$class<-"physiological-MTSV"

silv.funct.noint.dat<-full_join(extract_coefs(z.funct.drought.silvics.noint),extract_CIs(z.funct.drought.silvics.noint),by="trait")
colnames(silv.funct.noint.dat)<-c("trait","estimate","low","high")
silv.funct.noint.dat$class<-"functional-USFS"

silv.phys.noint.dat<-full_join(extract_coefs(z.phys.drought.silvics.noint),extract_CIs(z.phys.drought.silvics.noint),by="trait")
colnames(silv.phys.noint.dat)<-c("trait","estimate","low","high")
silv.phys.noint.dat$class<-"physiological-USFS"

mich.noint<-rbind(mich.phys.noint.dat,mich.funct.noint.dat)
USFS.noint<-rbind(silv.phys.noint.dat,silv.funct.noint.dat)

comps.noint<-rbind(mich.noint,USFS.noint)
comps.noint$category<-NA
comps.noint$category[which(comps.noint$class=="physiological-USFS")] <- "physiological"
comps.noint$category[which(comps.noint$class=="physiological-MTSV")] <- "physiological"
comps.noint$category[which(comps.noint$class=="functional-USFS")] <- "functional"
comps.noint$category[which(comps.noint$class=="functional-MTSV")] <- "functional"

comps.noint$data<-NA
comps.noint$data[which(comps.noint$class=="physiological-USFS")] <- "USFS"
comps.noint$data[which(comps.noint$class=="physiological-MTSV")] <- "MTSV"
comps.noint$data[which(comps.noint$class=="functional-USFS")] <- "USFS"
comps.noint$data[which(comps.noint$class=="functional-MTSV")] <- "MTSV"

###plotting==================================================
comps.noint$trait[which(comps.noint$trait=="pol_cent")] <- "pollination syndrome"
comps.noint$trait[which(comps.noint$trait=="flo_cent")] <- "flowering time"
comps.noint$trait[which(comps.noint$trait=="precip_cent")] <- "minimum precipitation"



pd=position_dodgev(height=0.4)
comps.noint.MTSV<-filter(comps.noint,data=="MTSV")
comps.noint.USFS<-filter(comps.noint,data=="USFS")
plotty4<-ggplot(comps.noint,aes(estimate,trait))+geom_point(size=2.5,aes(color=category,shape=data),position=pd)+geom_errorbarh(position=pd,width=0.4,aes(xmin=low,xmax=high,color=category,shape=data))+geom_vline(aes(xintercept=0))+theme_bw()+scale_color_manual(values=c("orchid4", "springgreen4"))

plotty4a<-ggplot(comps.noint.MTSV,aes(estimate,trait))+geom_point(size=2.5,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("MTSV")+xlim(-6.5,5)
plotty4b<-ggplot(comps.noint.USFS,aes(estimate,trait))+geom_point(size=2.5,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("USFS")+xlim(-6.5,5)


twopan.noint<-grid.arrange(plotty4a,plotty4b,nrow=1)
save.image(file="hystmodels.RData")


