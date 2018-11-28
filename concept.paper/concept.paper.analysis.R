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

#if you dont want to run the model: 
load("hystmodels.RData")
#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
drought.dat<-read.csv("..//Data/USDA_traitfor_MTSV.csv",header=TRUE)

###make a column for seed development time
mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time
###one more cleaninging tax
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)

###make the tree work 
mich.tree$node.label<-NULL

###RESPONSE VARIABLE KEY
#pro<- hysteranthy= before, and before with leaves 
#pro2<-hysteranthy= before, before/with and with
#pro3<- hysteranthy=before only
###phylo.D###### calculate phylo d.

set.seed(122)
###################################################################################
######phylo signals#########################################################
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloD <- phylo.d(d, binvar=pro) ###intermediate hysteranthy
PhyloD
plot(PhyloD)
Wind<-phylo.d(d, binvar=pol) #-0.48 more phylogenetically conserved that 
Wind
tol<-phylo.d(d, binvar=shade_bin)
tol
plot(Wind)
##functionalhysteranthy
PhyloPro2<-phylo.d(d,binvar=pro2)
PhyloPro2
plot(PhyloPro2)
#d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloPro3<-phylo.d(d,binvar=pro3)
PhyloPro3
plot(PhyloPro3)
###phlosignal for continuous trait
phylosig(mich.tree, mich.data$flo_time, method="lambda", test=TRUE, nsim=999,se=NULL)
phylosig(mich.tree, mich.data$dev.time, method="lambda", test=TRUE, nsim=999)
phylosig(mich.tree, mich.data$heigh_height, method="lambda", test=TRUE, nsim=999)
###############################################################################################

###Rescale predictors this makes it so you can compare binary to continous data
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))
mich.data$dev_time_cent<-(mich.data$dev.time-mean(mich.data$dev.time))/(2*sd(mich.data$dev.time))
mich.data$tol_cent<-(mich.data$shade_bin-mean(mich.data$shade_bin))/(2*sd(mich.data$shade_bin))



######### later we'll be running model with drought tolerance, so add drought tolerance data
mich.data<-left_join(mich.data,drought.dat,by="name")

######prune the tree for drought modeling
mich.data.wdrought<-filter(mich.data,!is.na(min._precip))
mich.data.wdrought$precip_cent<-(mich.data.wdrought$min._precip-mean(mich.data.wdrought$min._precip))/(2*sd(mich.data.wdrought$min._precip))
####prune tree to match reduced dataset
names.intree<-mich.tree$tip.label
namelist<-unique(mich.data.wdrought$name)
to.prune<-which(!names.intree%in%namelist)
mich.tree.droughtprune<-drop.tip(mich.tree,to.prune)
mytree.names<-mich.tree.droughtprune$tip.label
### the smaller new tree is called mich.tree.droughtprune

######### phyloglm requires species names to be in rownames
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")

###models:

#number of bootstraps from Wilcox, R. R. (2010). Fundamentals of modern statistical methods: Substantially improving power and accuracy. Springer.

######models with flowertime and polination and dev time, and complete species,
z.funct.seed<-phyloglm(pro2~pol_cent+flo_cent+pol_cent:flo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

######This code plots everything
bootestZ<-as.data.frame(z.funct.seed$coefficients)
bootconfZ<-as.data.frame(z.funct.seed$bootconfint95)
bootconfZ<-as.data.frame(t(bootconfZ))
bootestZ<-rownames_to_column(bootestZ, "trait")
bootconfZ<-rownames_to_column(bootconfZ, "trait")
bootmichZ<-full_join(bootconfZ,bootestZ, by="trait")
colnames(bootmichZ)<-c("trait","low","high","estimate")
bootmichZ<-dplyr::filter(bootmichZ, trait!="alpha")
bootmichZ<-dplyr::filter(bootmichZ, trait!="(Intercept)")
bootmichZ$class<-"functional-MTSV"

##############################################################
z.phys.seed<-phyloglm(pro3~pol_cent+flo_cent+pol_cent:flo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

###This code plots everything
bootestY<-as.data.frame(z.phys.seed$coefficients)
bootconfY<-as.data.frame(z.phys.seed$bootconfint95)
bootconfY<-as.data.frame(t(bootconfY))
bootestY<-rownames_to_column(bootestY, "trait")
bootconfY<-rownames_to_column(bootconfY, "trait")
bootmichY<-full_join(bootconfY,bootestY, by="trait")
colnames(bootmichY)<-c("trait","low","high","estimate")
bootmichY<-dplyr::filter(bootmichY, trait!="alpha")
bootmichY<-dplyr::filter(bootmichY, trait!="(Intercept)")
bootmichY$class<-"physiological-MTSV"

######### intermediate##############
z.inter.seed<-phyloglm(pro~pol_cent+flo_cent+pol_cent:flo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

bootestX<-as.data.frame(z.inter.seed$coefficients)
bootconfX<-as.data.frame(z.inter.seed$bootconfint95)
bootconfX<-as.data.frame(t(bootconfX))
bootestX<-rownames_to_column(bootestX, "trait")
bootconfX<-rownames_to_column(bootconfX, "trait")
bootmichX<-full_join(bootconfX,bootestX, by="trait")
colnames(bootmichX)<-c("trait","low","high","estimate")
bootmichX<-dplyr::filter(bootmichX, trait!="alpha")
bootmichX<-dplyr::filter(bootmichX, trait!="(Intercept)")
bootmichX$class<-"intermidiate-MTSV"
#######now put these all together
michcomp<-rbind(bootmichX,bootmichY,bootmichZ)

###and plot
pd=position_dodgev(height=0.25)
ggplot(michcomp,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()

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


###silv rescaling
silv.data$height_cent<-(silv.data$height-mean(silv.data$height))/(2*sd(silv.data$height))
silv.data$fruit_cent<-(silv.data$fruiting-mean(silv.data$fruiting))/(2*sd(silv.data$fruiting))
silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))
silv.data$pol_cent<-(silv.data$pol-mean(silv.data$pol))/(2*sd(silv.data$pol))
silv.data$tol_cent<-(silv.data$shade_bin-mean(silv.data$shade_bin))/(2*sd(silv.data$shade_bin))
silv.data$dev.time<-silv.data$fruiting-silv.data$flower_time
silv.data$dev_time_cent<-(silv.data$dev.time-mean(silv.data$dev.time))/(2*sd(silv.data$dev.time))

  
#######make a datasheet for dealing with slilvics drought
silv.data.wdrought<-filter(silv.data,!is.na(min._precip)) 
silv.data.wdrought$precip_cent<-(silv.data.wdrought$min._precip-mean(silv.data.wdrought$min._precip))/(2*sd(silv.data.wdrought$min._precip))    

####prune silvics tree to match reduced dataset
names.intree<-silv.tree$tip.label
namelist<-unique(silv.data.wdrought$name)
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
PhyloSilv3<-phylo.d(e,binvar=pro3)
PhyloSilv3
#######################

###prepare silvics data for modeling
silv.data<- silv.data %>% remove_rownames %>% column_to_rownames(var="name")

###models
z.funct.silv<-phyloglm(pro2~pol_cent*flo_cent,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

###prepare for plotting
bootestS<-as.data.frame(z.funct.silv$coefficients)
bootconfS<-as.data.frame(z.funct.silv$bootconfint95)
bootconfS<-as.data.frame(t(bootconfS))
bootestS<-rownames_to_column(bootestS, "trait")
bootconfS<-rownames_to_column(bootconfS, "trait")
bootsilvF<-full_join(bootconfS,bootestS, by="trait")
colnames(bootsilvF)<-c("trait","low","high","estimate")
bootsilvF<-dplyr::filter(bootsilvF, trait!="alpha")
bootsilvF<-dplyr::filter(bootsilvF, trait!="(Intercept)")
bootsilvF$class<-"functional-Silvics"

z.phys.silv<-phyloglm(pro3~pol_cent*flo_cent,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

bootestS1<-as.data.frame(z.phys.silv$coefficients)
bootconfS1<-as.data.frame(z.phys.silv$bootconfint95)
bootconfS1<-as.data.frame(t(bootconfS1))
bootestS1<-rownames_to_column(bootestS1, "trait")
bootconfS1<-rownames_to_column(bootconfS1, "trait")
bootsilvP<-full_join(bootconfS1,bootestS1, by="trait")
colnames(bootsilvP)<-c("trait","low","high","estimate")
bootsilvP<-dplyr::filter(bootsilvP, trait!="alpha")
bootsilvP<-dplyr::filter(bootsilvP, trait!="(Intercept)")
bootsilvP$class<-"physiological-Silvics"


z.inter.silv<-phyloglm(pro~pol_cent*flo_cent,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)

bootestSI<-as.data.frame(z.inter.silv$coefficients)
bootconfSI<-as.data.frame(z.inter.silv$bootconfint95)
bootconfSI<-as.data.frame(t(bootconfSI))
bootestSI<-rownames_to_column(bootestSI, "trait")
bootconfSI<-rownames_to_column(bootconfSI, "trait")
bootsilvI<-full_join(bootconfSI,bootestSI, by="trait")
colnames(bootsilvI)<-c("trait","low","high","estimate")
bootsilvI<-dplyr::filter(bootsilvI, trait!="alpha")
bootsilvI<-dplyr::filter(bootsilvI, trait!="(Intercept)")
bootsilvI$class<-"intermidiate-Silvics"

##combind into one data set
bootdatasilv<-rbind(bootsilvF,bootsilvP)
bootdatasilv<-rbind(bootdatasilv,bootsilvI)

##and plot
pd=position_dodgev(height=0.25)
ggplot(bootdatasilv,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()

###now combine mich and silvics
fullcomp<-rbind(michcomp,bootdatasilv)

##a little more formatting
fullcomp$category<-NA
fullcomp$category[which(fullcomp$class=="physiological-Silvics")] <- "physiological"
fullcomp$category[which(fullcomp$class=="physiological-MTSV")] <- "physiological"
fullcomp$category[which(fullcomp$class=="functional-Silvics")] <- "functional"
fullcomp$category[which(fullcomp$class=="functional-MTSV")] <- "functional"
fullcomp$category[which(fullcomp$class=="intermidiate-MTSV")] <- "intermediate"
fullcomp$category[which(fullcomp$class=="intermidiate-Silvics")] <- "intermediate"
fullcomp$data<-NA
fullcomp$data[which(fullcomp$class=="physiological-Silvics")] <- "USFS"
fullcomp$data[which(fullcomp$class=="physiological-MTSV")] <- "MTSV"
fullcomp$data[which(fullcomp$class=="functional-Silvics")] <- "USFS"
fullcomp$data[which(fullcomp$class=="functional-MTSV")] <- "MTSV"
fullcomp$data[which(fullcomp$class=="intermidiate-MTSV")] <- "MTSV"
fullcomp$data[which(fullcomp$class=="intermidiate-Silvics")] <- "USFS"



###and plot
pd=position_dodgev(height=0.6)
figure<-ggplot(fullcomp,aes(estimate,trait))+geom_point(size=2.5,aes(color=category,shape=data),position=pd,data=fullcomp)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category,group=class))+geom_vline(aes(xintercept=0))+xlim(-7.5,7)+theme_bw()
figure

###############################################################
##now add drought hypothesis... for mich
#################################################################
mich.data.wdrought<-  mich.data.wdrought %>% remove_rownames %>% column_to_rownames(var="name")

z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data.wdrought, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

z.phys.drought<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data.wdrought, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

z.inter.drought<-phyloglm(pro~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data.wdrought, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

bootestSI<-as.data.frame(z.inter.drought$coefficients)
bootconfSI<-as.data.frame(z.inter.drought$bootconfint95)
bootconfSI<-as.data.frame(t(bootconfSI))
bootestSI<-rownames_to_column(bootestSI, "trait")
bootconfSI<-rownames_to_column(bootconfSI, "trait")
bootdroughtI<-full_join(bootconfSI,bootestSI, by="trait")
colnames(bootdroughtI)<-c("trait","low","high","estimate")
bootdroughtI<-dplyr::filter(bootdroughtI, trait!="alpha")
bootdroughtI<-dplyr::filter(bootdroughtI, trait!="(Intercept)")
bootdroughtI$class<-"intermediate-MTSV"


bootestSF<-as.data.frame(z.funct.drought$coefficients)
bootconfSF<-as.data.frame(z.funct.drought$bootconfint95)
bootconfSF<-as.data.frame(t(bootconfSF))
bootestSF<-rownames_to_column(bootestSF, "trait")
bootconfSF<-rownames_to_column(bootconfSF, "trait")
bootdroughtF<-full_join(bootconfSF,bootestSF, by="trait")
colnames(bootdroughtF)<-c("trait","low","high","estimate")
bootdroughtF<-dplyr::filter(bootdroughtF, trait!="alpha")
bootdroughtF<-dplyr::filter(bootdroughtF, trait!="(Intercept)")
bootdroughtF$class<-"functional-MTSV"

bootestSP<-as.data.frame(z.phys.drought$coefficients)
bootconfSP<-as.data.frame(z.phys.drought$bootconfint95)
bootconfSP<-as.data.frame(t(bootconfSP))
bootestSP<-rownames_to_column(bootestSP, "trait")
bootconfSP<-rownames_to_column(bootconfSP, "trait")
bootdroughtP<-full_join(bootconfSP,bootestSP, by="trait")
colnames(bootdroughtP)<-c("trait","low","high","estimate")
bootdroughtP<-dplyr::filter(bootdroughtP, trait!="alpha")
bootdroughtP<-dplyr::filter(bootdroughtP, trait!="(Intercept)")
bootdroughtP$class<-"physiological-MTSV"

###combine these three catagories
bootdrought<-rbind(bootdroughtF,bootdroughtP)
bootdrought<-rbind(bootdrought,bootdroughtI)

#and plot
plotty<-ggplot(bootdrought,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()
plotty

#############################################################################################
#####silvics########  with drought 
##########################################################################################
silv.data.wdrought<- silv.data.wdrought %>% remove_rownames %>% column_to_rownames(var="name")

z.funct.drought.silvics<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data.wdrought, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

z.phys.drought.silvics<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data.wdrought, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                  start.beta=NULL, start.alpha=NULL,
                                  boot=599,full.matrix = TRUE)

z.inter.drought.silvics<-phyloglm(pro~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data.wdrought, silv.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                 start.beta=NULL, start.alpha=NULL,
                                 boot=599,full.matrix = TRUE)


bootestSI<-as.data.frame(z.inter.drought.silvics$coefficients)
bootconfSI<-as.data.frame(z.inter.drought.silvics$bootconfint95)
bootconfSI<-as.data.frame(t(bootconfSI))
bootestSI<-rownames_to_column(bootestSI, "trait")
bootconfSI<-rownames_to_column(bootconfSI, "trait")
bootdroughtII<-full_join(bootconfSI,bootestSI, by="trait")
colnames(bootdroughtII)<-c("trait","low","high","estimate")
bootdroughtII<-dplyr::filter(bootdroughtII, trait!="alpha")
bootdroughtII<-dplyr::filter(bootdroughtII, trait!="(Intercept)")
bootdroughtII$class<-"intermediate-USFS"


bootestSF<-as.data.frame(z.funct.drought.silvics$coefficients)
bootconfSF<-as.data.frame(z.funct.drought.silvics$bootconfint95)
bootconfSF<-as.data.frame(t(bootconfSF))
bootestSF<-rownames_to_column(bootestSF, "trait")
bootconfSF<-rownames_to_column(bootconfSF, "trait")
bootdroughtFF<-full_join(bootconfSF,bootestSF, by="trait")
colnames(bootdroughtFF)<-c("trait","low","high","estimate")
bootdroughtFF<-dplyr::filter(bootdroughtFF, trait!="alpha")
bootdroughtFF<-dplyr::filter(bootdroughtFF, trait!="(Intercept)")
bootdroughtFF$class<-"functional-USFS"

bootestSP<-as.data.frame(z.phys.drought.silvics$coefficients)
bootconfSP<-as.data.frame(z.phys.drought.silvics$bootconfint95)
bootconfSP<-as.data.frame(t(bootconfSP))
bootestSP<-rownames_to_column(bootestSP, "trait")
bootconfSP<-rownames_to_column(bootconfSP, "trait")
bootdroughtPP<-full_join(bootconfSP,bootestSP, by="trait")
colnames(bootdroughtPP)<-c("trait","low","high","estimate")
bootdroughtPP<-dplyr::filter(bootdroughtPP, trait!="alpha")
bootdroughtPP<-dplyr::filter(bootdroughtPP, trait!="(Intercept)")
bootdroughtPP$class<-"physiological-USFS"

bootdrought.USFS<-rbind(bootdroughtFF,bootdroughtPP)
bootdrought.USFS<-rbind(bootdrought.USFS,bootdroughtII)

plotty2<-ggplot(bootdrought.USFS,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()
plotty2
comps<-rbind(bootdrought,bootdrought.USFS)

comps$category<-NA
comps$category[which(comps$class=="physiological-USFS")] <- "physiological"
comps$category[which(comps$class=="physiological-MTSV")] <- "physiological"
comps$category[which(comps$class=="functional-USFS")] <- "functional"
comps$category[which(comps$class=="functional-MTSV")] <- "functional"
comps$category[which(comps$class=="intermediate-MTSV")] <- "intermediate"
comps$category[which(comps$class=="intermediate-USFS")] <- "intermediate"
comps$data<-NA
comps$data[which(comps$class=="physiological-USFS")] <- "USFS"
comps$data[which(comps$class=="physiological-MTSV")] <- "MTSV"
comps$data[which(comps$class=="functional-USFS")] <- "USFS"
comps$data[which(comps$class=="functional-MTSV")] <- "MTSV"
comps$data[which(comps$class=="intermediate-MTSV")] <- "MTSV"
comps$data[which(comps$class=="intermediate-USFS")] <- "USFS"


plotty3<-ggplot(comps,aes(estimate,trait))+geom_point(size=2.5,aes(color=category,shape=data),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category,shape=data))+geom_vline(aes(xintercept=0))+theme_bw()
plotty3

save.image(file="hystmodels.RData")


####try a order catagorical model
unique(mich.data.wdrought$Phen.sequence)
mich.data.wdrought$hyst_cat<-NA
mich.data.wdrought$hyst_cat[which(mich.data.wdrought$Phen.sequence=="ser")]<-0
mich.data.wdrought$hyst_cat[which(mich.data.wdrought$Phen.sequence=="hyst")]<-0
mich.data.wdrought$hyst_cat[which(mich.data.wdrought$Phen.sequence=="syn/ser")]<-1
mich.data.wdrought$hyst_cat[which(mich.data.wdrought$Phen.sequence=="syn")]<-2
mich.data.wdrought$hyst_cat[which(mich.data.wdrought$Phen.sequence=="pro/syn")]<-3
mich.data.wdrought$hyst_cat[which(mich.data.wdrought$Phen.sequence=="pro")]<-4
unique(mich.data.wdrought$hyst_cat)



mod.cat<-polr(factor(hyst_cat)~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent, mich.data.wdrought)
summary(mod.cat)
exp(coef(mod.cat))


z.funct.cat<--phyloglm(factor(hyst_cat)~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data.wdrought, mich.tree.droughtprune, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)
