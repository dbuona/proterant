rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant")
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
library(lme4)
source("Scripts/continuous_mod_prep.R")
###individuals with more than 10 observations
checker1<-d %>% group_by(tree.id) %>% summarise(non_na_count = sum(!is.na(offset.funct)))
checker1<-filter(checker1, non_na_count>=10)
table(checker1$tree.id)
use.id<-checker1$tree.id
d.plus<-filter(d, tree.id %in% use.id)

##species to use
###ACRU, QURU, ACPE, most complete observations, each hysteranthy class

###species with more than 45 total boservations
use.sp<-c("ACPE","ACRU","QURU")
#"KALA","BEAL","AMSP","VACO

d.plus<-filter(d.plus, species %in% use.sp)

jpeg("../figure/individual_HF_var.jpg")
ggplot(d.plus)+stat_summary(fun.data = "mean_cl_boot",geom="errorbar",aes(tree.id, offset.phys,color=species),linetype="solid")+stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",aes(tree.id,offset.funct,color=species),linetype="dashed")+stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",aes(tree.id,offset.inter,color=species),linetype="dotdash")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+ylab("FLS offset")+xlab("tree ID")
dev.off()

m<-lmer(offset.inter~tree.id+(1|species)+(1|year),data=d.plus)
n<-lmer(offset.phys~tree.id+(1|species)+(1|year),data=d.plus)
o<-lmer(offset.funct~tree.id+(1|species)+(1|year),data=d.plus)
Anova(m)
Anova(n)
Anova(o)
####This indicates that idividuals vary in the phys and inter, but bot
#### interannual

jpeg("../figure/interannual_HF_var.jpg")
ggplot(d.plus)+geom_line(aes(year,offset.funct, group=tree.id, color=species),linetype="dashed")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+geom_line(aes(year,offset.phys, group=tree.id, color=species),linetype="solid")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+geom_line(aes(year,offset.inter, group=tree.id, color=species),linetype="dotdash")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+ylab("FLS offset")
dev.off()

mm<-lmer(offset.inter~year+(1|tree.id),data=d.plus)
nn<-lmer(offset.phys~year+(1|tree.id),data=d.plus)
oo<-lmer(offset.funct~year+(1|tree.id),data=d.plus)
summary(mm)
Anova(mm)


  ###This is the continuous model not using michigan tree so we can address other hypothesese







setdiff(d$species, traits$species)
intersect(d$species, traits$species)

spfordata<-traits$species

d<-filter(d,species %in% c(spfordata))

dater<-left_join(d,traits, by="species")

intersect(traits$species, dater$species)
unique(dater$species)

#flower_ave<-dater %>% group_by(name) %>% summarise(average.flower.time=mean(fopn.jd,na.rm=TRUE))
flower_ave<-dplyr::select(dater,tree.id,species,name, year,fopn.jd)

flower_ave$origin <- as.Date(paste0(flower_ave$year, "-01-01"),tz = "UTC") 
flower_ave$flo_date<-as.Date(flower_ave$fopn.jd, origin = flower_ave$origin, tz = "UTC") 
flower_ave<-separate(flower_ave, flo_date, c("year", "month", "dia"))
flower_ave$month<-as.numeric(flower_ave$month)

flo_month<-flower_ave %>% group_by(name) %>% summarise(flowering_month=mean(month,na.rm=TRUE))


dater<-left_join(dater,flo_month, by="name")
table(dater$name)

dater$cent_flo_month<-(dater$flowering_month-mean(dater$flowering_month,na.rm=TRUE))/(2*sd(dater$flowering_month,na.rm=TRUE))


library(brms)
library(MCMCglmm)


inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###bayesian and continuous-- main model###############
###all measurements of flowering time swamp all other predictors
modelcont.funct <- brm(offset.funct~ cent_pol+seed_cent+cent_minP+cent_minT+cent_flo_month +(1|name), data = dater, 
                 family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelcont.phys <- brm(offset.phys~ cent_pol+seed_cent+cent_minP+cent_minT+cent_flo_month+(1|name), data = dater, 
                 family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelcont.inter <- brm(offset.inter~ cent_pol+seed_cent+cent_minP+cent_minT+cent_flo_month+(1|name), data = dater, 
                      family = gaussian(), cov_ranef = list(name= A),iter=3000) 

A<-as.data.frame(tidy(modelcont,robust = TRUE))
A<-A %>% "["(.,2:6,)
A$class<-"functional"

B<-as.data.frame(tidy(modelcont.inter,robust = TRUE))
B<-B %>% "["(.,2:6,)
B$class<-"intermediate"


C<-as.data.frame(tidy(modelcont.phys,robust = TRUE))
C<-C %>% "["(.,2:6,)
C$class<-"physiological"
D<-rbind(A,B,C)

library(ggstance)
pd=position_dodgev(height=0.3)
ggplot(D,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="black")
######## Do models again without flowering time
colnames(dater)
modelcont.funct.1 <- brm(offset.funct~ cent_pol+seed_cent+cent_minP+cent_minT+cent_flo_view+(1|name), data = dater, 
                 family = gaussian(), cov_ranef = list(name= A),iter=3000) 

summary(modelcont.funct.1)

modelcont.phys.1 <- brm(offset.phys~ cent_pol+(1|name), data = dater, 
family = gaussian(), cov_ranef = list(name= A),iter=3000) 

summary(modelcont.phys.1)

modelcont.inter.1 <- brm(offset.inter~ cent_pol+seed_cent+cent_minP+cent_minT+cent_frost+cent_roots +(1|name), data = dater, 
                       family = gaussian(), cov_ranef = list(name= A),iter=3000) 

coef(modelcont.funct.1)
A<-as.data.frame(tidy(modelcont.funct.1,robust = TRUE))
A<-A %>% "["(.,2:7,)
A$class<-"functional"

B<-as.data.frame(tidy(modelcont.inter.1,robust = TRUE))
B<-B %>% "["(.,2:7,)
B$class<-"intermediate"


C<-as.data.frame(tidy(modelcont.phys.1,robust = TRUE))
C<-C %>% "["(.,2:7,)
C$class<-"physiological"
D<-rbind(A,B,C)

library(ggstance)
pd=position_dodgev(height=0.3)
ggplot(D,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="black")
