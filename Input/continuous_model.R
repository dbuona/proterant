###run a model using harvard forest species combine with michigan trees data
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

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time
###one more cleaninging tax
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)

###Rescale predictors
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$dev_time_cent<-(mich.data$dev.time-mean(mich.data$dev.time))/(2*sd(mich.data$dev.time))
mich.data$tol_cent<-(mich.data$shade_bin-mean(mich.data$shade_bin))/(2*sd(mich.data$shade_bin))

##read in harvard forest
d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
unique(d$species)
### offset functional
d$offset<-d$l75.jd-d$fopn.jd
###0ffset physiological
d$offset2<-d$bb.jd-d$fbb.jd
###intermediate
d$offset3<-d$bb.jd-d$fopn.jd

########combine datasets
d$name[d$species=="ACPE"]<-"Acer_pensylvanicum"
d$name[d$species=="ACRU"]<-"Acer_rubrum"
d$name[d$species=="ACSA"]<-"Acer_saccharum"
d$name[d$species=="AMSP"]<-"Amelanchier_arborea"
d$name[d$species=="BEAL"]<-"Betula_alleghaniensis"
d$name[d$species=="BEPA"]<-"Betula_papyrifera"
d$name[d$species=="FAGR"]<-"Fagus_grandifolia"
d$name[d$species=="FRAM"]<-"Fraxinus_americana"
d$name[d$species=="HAVI"]<-"Hamamelis_virginiana"
d$name[d$species=="ILVE"]<-"Ilex_verticillata"
d$name[d$species=="KAAN"]<-"Kalmia_angustifolia"
d$name[d$species=="NEMU"]<-"Ilex_mucronata"
d$name[d$species=="NYSY"]<-"Nyssa_sylvatica"
d$name[d$species=="POTR"]<-"Populus_tremuloides"
d$name[d$species=="PRSE"]<-"Prunus_serotina"
d$name[d$species=="QURU"]<-"Quercus_rubra"
d$name[d$species=="QUVE"]<-"Quercus_velutina"
d$name[d$species=="VACO"]<-"Vaccinium_corymbosum"
d$name[d$species=="VICA"]<-"Viburnum_cassinoides"
d$name[d$species=="SAPU"]<-"Sambucus_racemosa"
d$name[d$species=="COAL"]<-"Cornus_alternifolia"

cont<-left_join(d,mich.data,by="name")
cont<-filter(cont,!is.na(offset))
cont<-filter(cont,!is.na(name))
unique(cont$name)
setdiff(cont$name,d$name)

names.intree<-mich.tree$tip.label

# list of my species myspecies
namelist<-unique(d$name)
####clean
setdiff(d$name,namelist)

to.prune<-which(!names.intree%in%namelist)
HF.tree<-drop.tip(mich.tree,to.prune)

sum <-cont %>% group_by(name) %>% summarise(avg.offset = mean(offset3,na.rm=TRUE))

cont.dat<-left_join(sum,mich.data)

###format the data in the same order as the tree
mytree.names<-HF.tree$tip.label

final.df<-cont.dat[match(mytree.names, cont.dat$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names
HF.tree$node.label<-NULL


final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")

plot(HF.tree)
is.ultrametric(HF.tree)

HF.mod<-lm(avg.offset~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,final.df)
summary(HF.mod)
fit <- phylolm(avg.offset~pol+heigh_height+flo_time+dev.time+shade_bin,data=final.df,phy=HF.tree,model="BM",measurement_error=TRUE,boot=599)

summary(fit)

bootest2<-as.data.frame(fit$coefficients)
bootconf2<-as.data.frame(fit$bootconfint95)
bootconf2<-as.data.frame(t(bootconf2))

bootest2<-rownames_to_column(bootest2, "trait")
bootconf2<-rownames_to_column(bootconf2, "trait")
bootmich2<-full_join(bootconf2,bootest2, by="trait")
colnames(bootmich2)<-c("trait","low","high","estimate")
bootmich2<-dplyr::filter(bootmich2, trait!="sigma2")
bootmich2<-dplyr::filter(bootmich2, trait!="sigma2_error")
bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")

library(ggthemes)
functplot1<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-20,60)+theme_base()
functplot1

###scaled model
fit2 <- phylolm(avg.offset~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,phy=HF.tree,model="BM",measurement_error=TRUE,boot=599)
summary(fit2)
bootest2<-as.data.frame(fit2$coefficients)
bootconf2<-as.data.frame(fit2$bootconfint95)
bootconf2<-as.data.frame(t(bootconf2))

bootest2<-rownames_to_column(bootest2, "trait")
bootconf2<-rownames_to_column(bootconf2, "trait")
bootmich2<-full_join(bootconf2,bootest2, by="trait")
colnames(bootmich2)<-c("trait","low","high","estimate")
bootmich2<-dplyr::filter(bootmich2, trait!="sigma2")
bootmich2<-dplyr::filter(bootmich2, trait!="sigma2_error")
bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")

functplot2<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-60,60)+theme_base()
functplot2

summary(fit)
unique(cont.dat$name)

###brms

sum2 <-cont %>% group_by(name,tree.id) %>% summarise(avg.offset = mean(offset3,na.rm=TRUE))

cont.dat2<-left_join(sum2,mich.data)

###format the data in the same order as the tree
sum2 <-cont %>% group_by(name,tree.id) %>% summarise(avg.offset = mean(offset3,na.rm=TRUE))

cont.dat2<-left_join(sum2,mich.data)

final.df2<-cont.dat2[match(mytree.names, cont.dat2$name),]
namelist<-final.df2$name
namelist==mytree.names
final.df2$name== mytree.names
HF.tree$node.label<-NULL

library("MCMCglmm")
library(brms)
final.dff<-rownames_to_column(cont.dat, "name")
inv.phylo <- MCMCglmm::inverseA(HF.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###bayesian and continuous-- main model###############
modelcont.brms <- brm(avg.offset~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent+(1|name), data = cont.dat, 
                 family = gaussian(), cov_ranef = list(name= A),warmup=2000,iter=3000) 

summary(modelcont.brms)
summary(fit2)

prior = c(prior(normal(0, 5), "b"),
          prior(normal(0, 5), "Intercept"),
          prior(student_t(3, 0, 5), "sd"))


