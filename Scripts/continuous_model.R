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
library(broom)
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

#load(file="contmodels.RData")
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

###make column for mich trees rownames
mich.data<-rownames_to_column(mich.data, "name")

cont<-left_join(d,mich.data,by="name")
cont<-filter(cont,!is.na(offset))
cont<-filter(cont,!is.na(name))
unique(cont$name)
setdiff(cont$name,d$name)
###Binary
cont$hyst<-ifelse(cont$offset>0,1,0)
cont$hyst2<-ifelse(cont$offset2>0,1,0)
cont$hyst3<-ifelse(cont$offset2>0,1,0)


names.intree<-mich.tree$tip.label

# list of my species myspecies
namelist<-unique(d$name)
####clean
setdiff(d$name,namelist)

to.prune<-which(!names.intree%in%namelist)
HF.tree<-drop.tip(mich.tree,to.prune)

####

inv.phylo <- MCMCglmm::inverseA(HF.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###bayesian and continuous-- main model###############
library(brms)
modelcont <- brm(offset~ pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent +(1|name), data = cont, 
                 family = gaussian, cov_ranef = list(name= A) ,iter=5000)                
summary(modelcont)

modelcont2 <- brm(offset2~ pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent +(1|name), data = cont, 
                 family = gaussian, cov_ranef = list(name= A) ,iter=5000)                
summary(modelcont2)

modelcont3 <- brm(offset3~ pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent +(1|name), data = cont, 
                  family = gaussian, cov_ranef = list(name= A) ,iter=5000)                
summary(modelcont3)

modelbin<-brm(hyst~ pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent +(1|name), data = cont, 
              family = bernoulli(link="logit"), cov_ranef = list(name= A) ,iter=5000) 

modelbin2<-brm(hyst2~ pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent +(1|name), data = cont, 
              family = bernoulli(link="logit"), cov_ranef = list(name= A) ,iter=5000) 

modelbin3<-brm(hyst3~ pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent +(1|name), data = cont, 
              family = bernoulli(link="logit"), cov_ranef = list(name= A) ,iter=5000) 

A<-as.data.frame(tidy(modelcont,robust = TRUE))
A<-A %>% "["(.,2:6,)
A$class<-"functional"
  
B<-as.data.frame(tidy(modelcont2,robust = TRUE))
B<-B %>% "["(.,2:6,)
B$class<-"intermediate"


C<-as.data.frame(tidy(modelcont3,robust = TRUE))
C<-C %>% "["(.,2:6,)
C$class<-"physiological"

AA<-as.data.frame(tidy(modelbin,robust = TRUE))
AA<-AA %>% "["(.,2:6,)
AA$class<-"functional"

BB<-as.data.frame(tidy(modelbin2,robust = TRUE))
BB<-BB %>% "["(.,2:6,)
BB$class<-"intermediate"

CC<-as.data.frame(tidy(modelbin3,robust = TRUE))
CC<-CC %>% "["(.,2:6,)
CC$class<-"physiological"

D<-rbind(A,B,C)
DD<-rbind(AA,BB,CC)
library(ggstance)
pd=position_dodgev(height=0.3)
ggplot(D,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="blue")

ggplot(DD,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="blue")






######go here to run the model in phylglm
sum1 <-cont %>% group_by(name) %>% summarise(avg.offset = mean(offset3,na.rm=TRUE)) 
sum2 <-cont %>% group_by(name) %>% summarise(avg.offset.phys = mean(offset2,na.rm=TRUE))
sum3 <-cont %>% group_by(name) %>% summarise(avg.offset.func = mean(offset,na.rm=TRUE))

sum<-left_join(sum1,sum2,by="name")
sum<-left_join(sum,sum3,by="name")
cont.dat<-left_join(sum,mich.data,by="name")

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

######un scaled#######################################################
#fit <- phylolm(avg.offset~pol+heigh_height+flo_time+dev.time+shade_bin,data=final.df,phy=HF.tree,model="BM",measurement_error=TRUE,boot=599)
#summary(fit)
#bootest2<-as.data.frame(fit$coefficients)
#bootconf2<-as.data.frame(fit$bootconfint95)
#bootconf2<-as.data.frame(t(bootconf2))

#bootest2<-rownames_to_column(bootest2, "trait")
#bootconf2<-rownames_to_column(bootconf2, "trait")
#bootmich2<-full_join(bootconf2,bootest2, by="trait")
#colnames(bootmich2)<-c("trait","low","high","estimate")
#bootmich2<-dplyr::filter(bootmich2, trait!="sigma2")
#bootmich2<-dplyr::filter(bootmich2, trait!="sigma2_error")
#bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")

#library(ggthemes)
#functplot1<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-20,60)+theme_base()
#functplot1

###scaled model
fit2 <- phylolm(avg.offset~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,phy=HF.tree,model="BM",measurement_error=TRUE,boot=599)
summary(fit2)
bootest<-as.data.frame(fit2$coefficients)
bootconf<-as.data.frame(fit2$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="sigma2")
bootmich<-dplyr::filter(bootmich, trait!="sigma2_error")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")

bootmich$trait[bootmich$trait=="tol_cent"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol_cent"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="dev_time_cent"]<-"seed development"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

###take reciprocal of flowering timg
bootmich$low<-ifelse(bootmich$trait==c("flower timing"),-(bootmich$low),bootmich$low)
bootmich$high<-ifelse(bootmich$trait==c("flower timing"),-(bootmich$high),bootmich$high)
bootmich$estimate<-ifelse(bootmich$trait==c("flower timing"),-(bootmich$estimate),bootmich$estimate)

bootmich$trait[bootmich$trait=="flower timing"]<-"earlier flowering"
bootmich$trait[bootmich$trait=="seed development"]<-"seed development period"
bootmich$data<-"MTSV"
bootmich$category<-"intermediate"

contplot<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5,aes(shape=data,color=category))+geom_segment(aes(y=trait,yend=trait,x=low,xend=high,color=category))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0),color="black")+xlim(-58,60)+theme_bw()+annotate("text", x = 52, y = 5.5, label = "Flowers first",fontface =2)+annotate("text", x = -51, y = 5.5, label = "Leaves first",fontface =2)
contplot
library("sme")
AIC(fit2)
AIC(fit2.func)
AIC(fit2.physi)
AIC(funct.HF)
AIC(phys.HF)
AIC(inter.HF)
#####other form of hysteranthy: physiological
fit2.func <- phylolm(avg.offset.func~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,phy=HF.tree,model="BM",measurement_error=TRUE,boot=599)
summary(fit2.func)
bootest2<-as.data.frame(fit2.func$coefficients)
bootconf2<-as.data.frame(fit2.func$bootconfint95)
bootconf2<-as.data.frame(t(bootconf2))

bootest2<-rownames_to_column(bootest2, "trait")
bootconf2<-rownames_to_column(bootconf2, "trait")
bootmich2<-full_join(bootconf2,bootest2, by="trait")
colnames(bootmich2)<-c("trait","low","high","estimate")
bootmich2<-dplyr::filter(bootmich2, trait!="sigma2")
bootmich2<-dplyr::filter(bootmich2, trait!="sigma2_error")
bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")

bootmich2$trait[bootmich2$trait=="tol_cent"]<-"shade tolerance"
bootmich2$trait[bootmich2$trait=="pol_cent"]<-"pollination syndrome"
bootmich2$trait[bootmich2$trait=="height_cent"]<-"max height"
bootmich2$trait[bootmich2$trait=="dev_time_cent"]<-"seed development"
bootmich2$trait[bootmich2$trait=="flo_cent"]<-"flower timing"

###take reciprocal of flowering timg
bootmich2$low<-ifelse(bootmich2$trait==c("flower timing"),-(bootmich2$low),bootmich2$low)
bootmich2$high<-ifelse(bootmich2$trait==c("flower timing"),-(bootmich2$high),bootmich2$high)
bootmich2$estimate<-ifelse(bootmich2$trait==c("flower timing"),-(bootmich2$estimate),bootmich2$estimate)

bootmich2$trait[bootmich2$trait=="flower timing"]<-"earlier flowering"
bootmich2$trait[bootmich2$trait=="seed development"]<-"seed development period"
bootmich2$data<-"MTSV"
bootmich2$category<-"functional"
contplot.func<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5,aes(color=category,shape=data))+geom_segment(aes(y=trait,yend=trait,x=low,xend=high,color=category))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0),color="black")+xlim(-58,60)+theme_bw()+annotate("text", x = 52, y = 5.5, label = "Flowers first",fontface =2)+annotate("text", x = -51, y = 5.5, label = "Leaves first",fontface =2)
contplot.func

#####physiological
fit2.physi <- phylolm(avg.offset.phys~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,phy=HF.tree,model="BM",measurement_error=TRUE,boot=599)
summary(fit2.physi)

bootestS1<-as.data.frame(fit2.physi$coefficients)
bootconfS1<-as.data.frame(fit2.physi$bootconfint95)

bootconfS1<-as.data.frame(t(bootconfS1))


bootestS1<-rownames_to_column(bootestS1, "trait")
bootconfS1<-rownames_to_column(bootconfS1, "trait")
bootsilvP<-full_join(bootconfS1,bootestS1, by="trait")
colnames(bootsilvP)<-c("trait","low","high","estimate")
bootsilvP<-dplyr::filter(bootsilvP, trait!="sigma2")
bootsilvP<-dplyr::filter(bootsilvP, trait!="sigma2_error")
bootsilvP<-dplyr::filter(bootsilvP, trait!="(Intercept)")
###names
bootsilvP$trait[bootsilvP$trait=="tol_cent"]<-"shade tolerance"
bootsilvP$trait[bootsilvP$trait=="pol_cent"]<-"pollination syndrome"
bootsilvP$trait[bootsilvP$trait=="height_cent"]<-"max height"
bootsilvP$trait[bootsilvP$trait=="dev_time_cent"]<-"seed development"
bootsilvP$trait[bootsilvP$trait=="flo_cent"]<-"flower timing"
bootsilvP$category<-"physiological"
bootsilvP$data<-"MTSV"

bootsilvP$low<-ifelse(bootsilvP$trait==c("flower timing"),-(bootsilvP$low),bootsilvP$low)
bootsilvP$high<-ifelse(bootsilvP$trait==c("flower timing"),-(bootsilvP$high),bootsilvP$high)
bootsilvP$estimate<-ifelse(bootsilvP$trait==c("flower timing"),-(bootsilvP$estimate),bootsilvP$estimate)

bootsilvP$trait[bootsilvP$trait=="flower timing"]<-"earlier flowering"
bootsilvP$trait[bootsilvP$trait=="seed development"]<-"seed development period"

contplot.physi<-ggplot(bootsilvP,aes(estimate,trait))+geom_point(size=2.5,aes(color=category,shape=data))+geom_segment(aes(y=trait,yend=trait,x=low,xend=high,color=category))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0),color="black")+xlim(-58,60)+theme_bw()+annotate("text", x = 52, y = 5.5, label = "Flowers first",fontface =2)+annotate("text", x = -51, y = 5.5, label = "Leaves first",fontface =2)
contplot.physi

full.cont<-rbind(bootmich,bootmich2,bootsilvP)


pd=position_dodgev(height=0.4)
contplot.full<-ggplot(full.cont,aes(estimate,trait))+geom_point(size=2.5,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,height=0,aes(xmin=low,xmax=high,color=category))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0),color="black")+xlim(-58,60)+theme_bw()+annotate("text", x = 52, y = 5.5, label = "Flowers first",fontface =2)+annotate("text", x = -51, y = 5.5, label = "Leaves first",fontface =2)
contplot.full
##########now do harvard forest as binary
final.df$hyst.inter<-ifelse(final.df$avg.offset>0,1,0)
final.df$hyst.func<-ifelse(final.df$avg.offset.func>0,1,0)
final.df$hyst.phys<-ifelse(final.df$avg.offset.phys>0,1,0)

######## binary models
funct.HF<-phyloglm(hyst.func~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,final.df, HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                      start.beta=NULL, start.alpha=NULL,
                      boot=10,full.matrix = TRUE)
summary(funct.HF)
bootestA<-as.data.frame(funct.HF$coefficients)
bootconfA<-as.data.frame(funct.HF$bootconfint95)
bootconfA<-as.data.frame(t(bootconfA))

bootestA<-rownames_to_column(bootestA, "trait")
bootconfA<-rownames_to_column(bootconfA, "trait")
bootmichA<-full_join(bootconfA,bootestA, by="trait")
colnames(bootmichA)<-c("trait","low","high","estimate")
bootmichA<-dplyr::filter(bootmichA, trait!="sigma2")
bootmichA<-dplyr::filter(bootmichA, trait!="sigma2_error")
bootmichA<-dplyr::filter(bootmichA, trait!="(Intercept)")

bootmichA$trait[bootmichA$trait=="tol_cent"]<-"shade tolerance"
bootmichA$trait[bootmichA$trait=="pol_cent"]<-"pollination syndrome"
bootmichA$trait[bootmichA$trait=="height_cent"]<-"max height"
bootmichA$trait[bootmichA$trait=="dev_time_cent"]<-"seed development"
bootmichA$trait[bootmichA$trait=="flo_cent"]<-"flower timing"
bootmichA$category<-"functional"
phys.HF<-phyloglm(hyst.phys~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,final.df, HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                   start.beta=NULL, start.alpha=NULL,
                   boot=10,full.matrix = TRUE)
summary(phys.HF)
bootestB<-as.data.frame(phys.HF$coefficients)
bootconfB<-as.data.frame(phys.HF$bootconfint95)
bootconfB<-as.data.frame(t(bootconfB))

bootestB<-rownames_to_column(bootestB, "trait")
bootconfB<-rownames_to_column(bootconfB, "trait")
bootmichB<-full_join(bootconfB,bootestB, by="trait")
colnames(bootmichB)<-c("trait","low","high","estimate")
bootmichB<-dplyr::filter(bootmichB, trait!="sigma2")
bootmichB<-dplyr::filter(bootmichB, trait!="sigma2_error")
bootmichB<-dplyr::filter(bootmichB, trait!="(Intercept)")

bootmichB$trait[bootmichB$trait=="tol_cent"]<-"shade tolerance"
bootmichB$trait[bootmichB$trait=="pol_cent"]<-"pollination syndrome"
bootmichB$trait[bootmichB$trait=="height_cent"]<-"max height"
bootmichB$trait[bootmichB$trait=="dev_time_cent"]<-"seed development"
bootmichB$trait[bootmichB$trait=="flo_cent"]<-"flower timing"
bootmichB$category<-"physiological"
inter.HF<-phyloglm(hyst.inter~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,final.df, HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                  start.beta=NULL, start.alpha=NULL,
                  boot=10,full.matrix = TRUE) ###THis worked
summary(inter.HF)
bootestC<-as.data.frame(inter.HF$coefficients)
bootconfC<-as.data.frame(inter.HF$bootconfint95)
bootconfC<-as.data.frame(t(bootconfC))

bootestC<-rownames_to_column(bootestC, "trait")
bootconfC<-rownames_to_column(bootconfC, "trait")
bootmichC<-full_join(bootconfC,bootestC, by="trait")
colnames(bootmichC)<-c("trait","low","high","estimate")
bootmichC<-dplyr::filter(bootmichC, trait!="sigma2")
bootmichC<-dplyr::filter(bootmichC, trait!="sigma2_error")
bootmichC<-dplyr::filter(bootmichC, trait!="(Intercept)")

bootmichC$trait[bootmichC$trait=="tol_cent"]<-"shade tolerance"
bootmichC$trait[bootmichC$trait=="pol_cent"]<-"pollination syndrome"
bootmichC$trait[bootmichC$trait=="height_cent"]<-"max height"
bootmichC$trait[bootmichC$trait=="dev_time_cent"]<-"seed development"
bootmichC$trait[bootmichC$trait=="flo_cent"]<-"flower timing"
bootmichC$category<-"intermediate"
HF.bin<-rbind(bootmichA,bootmichB,bootmichC)
HF.bin<-filter(HF.bin, trait!="alpha")

HF.bin$low<-ifelse(HF.bin$trait==c("flower timing"),-(HF.bin$low),HF.bin$low)
HF.bin$high<-ifelse(HF.bin$trait==c("flower timing"),-(HF.bin$high),HF.bin$high)
HF.bin$estimate<-ifelse(HF.bin$trait==c("flower timing"),-(HF.bin$estimate),HF.bin$estimate)

HF.bin$trait[HF.bin$trait=="flower timing"]<-"earlier flowering"
HF.bin$trait[HF.bin$trait=="seed development"]<-"seed development period"

pd=position_dodgev(height=0.4)
binplot.full<-ggplot(HF.bin,aes(estimate,trait))+geom_point(size=2.5,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,height=0,aes(xmin=low,xmax=high,color=category))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0),color="black")+xlim(-5,5.5)+theme_bw()+annotate("text", x = 5, y = 5.5, label = "Flowers first",fontface =2)+annotate("text", x = -4.2, y = 5.5, label = "Leaves first",fontface =2)
binplot.full

library(gridExtra)
product.plot<-grid.arrange(binplot.full,contplot.full, ncol=2)

save.image(file="contmodels.RData")
?phyloglm()
phys<-binaryPGLMM(hyst.phys~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,HF.tree) 
inter<-binaryPGLMM(hyst.inter~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,HF.tree)
func<-binaryPGLMM(hyst.func~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,data=final.df,HF.tree)
