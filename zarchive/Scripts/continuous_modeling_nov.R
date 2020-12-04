rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant")
library("phytools")
library(broom)
library(dplyr)
library(tidyr)
library(boot)
library("ggplot2")
library(car)
library(lme4)
library(brms)
library(ggstance)
library(MCMCglmm)


source("Scripts/continuous_mod_prep.R")
load("continuous.mods.RData")

##species to use
###ACRU, QURU, ACPE, most complete observations, each hysteranthy class

use.sp<-c("ACPE","ACRU","QURU")

d.plus<-filter(d, species %in% use.sp)

#jpeg("../figure/individual_HF_var.jpg")
#ggplot(d.plus)+stat_summary(fun.data = "mean_cl_boot",geom="errorbar",aes(tree.id, offset.phys,color=species),linetype="solid")+stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",aes(tree.id,offset.funct,color=species),linetype="dashed")+stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",aes(tree.id,offset.inter,color=species),linetype="dotdash")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+ylab("FLS offset")+xlab("tree ID")
#dev.off()

m<-lmer(offset.inter~tree.id+(1|species)+(1|year),data=d.plus)
n<-lmer(offset.phys~tree.id+(1|species)+(1|year),data=d.plus)
o<-lmer(offset.funct~tree.id+(1|species)+(1|year),data=d.plus)
Anova(m)
Anova(n)
Anova(o)
####This indicates that idividuals vary in the phys and inter, but bot
#### interannual

#jpeg("../figure/interannual_HF_var.jpg")
#ggplot(d.plus)+geom_line(aes(year,offset.funct, group=tree.id, color=species),linetype="dashed")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+geom_line(aes(year,offset.phys, group=tree.id, color=species),linetype="solid")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+geom_line(aes(year,offset.inter, group=tree.id, color=species),linetype="dotdash")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+ylab("FLS offset")
#dev.off()

mm<-lmer(offset.inter~year+(1|tree.id),data=d.plus)
nn<-lmer(offset.phys~year+(1|tree.id),data=d.plus)
oo<-lmer(offset.funct~year+(1|tree.id),data=d.plus)
summary(mm)
Anova(mm)


  ###This is the continuous model not using michigan tree so we can address other hypothesese


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
####recenter predictors based on new data
dater$cent_flo_month<-(dater$flowering_month-mean(dater$flowering_month,na.rm=TRUE))/(2*sd(dater$flowering_month,na.rm=TRUE))
dater$cent_pol<-(dater$pol-mean(traits$pol,na.rm=TRUE))/(2*sd(dater$pol,na.rm=TRUE))
dater$cent_minP<-(dater$min_precip-mean(dater$min_precip))/(2*sd(dater$min_precip))
dater$cent_floday<-(dater$fopn.jd-mean(dater$fopn.jd,na.rm=TRUE))/(2*sd(dater$fopn.jd,na.rm=TRUE))



####recenter predictors based on new data

library(brms)


inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###bayesian and continuous-- no intx###############
###all measurements of flowering time swamp all other predictors
#modelcont.funct <- brm(offset.funct~cent_pol+cent_minP+cent_flo_month +(1|name), data = dater, 
 #                family = gaussian(), cov_ranef = list(name= A),iter=3000) 

#modelcont.phys <- brm(offset.phys~ cent_pol+cent_minP+cent_flo_month+(1|name), data = dater, 
 #                family = gaussian(), cov_ranef = list(name= A),iter=3000) 

#modelcont.inter <- brm(offset.inter~ cent_pol+cent_minP+cent_flo_month+(1|name), data = dater, 
  #                    family = gaussian(), cov_ranef = list(name= A),iter=3000) 

#Aa<-as.data.frame(tidy(modelcont.funct,robust = TRUE))
#Aa<-Aa %>% "["(.,2:5,)
#Aa$class<-"functional"

#Bb<-as.data.frame(tidy(modelcont.inter,robust = TRUE))
#Bb<-Bb %>% "["(.,2:5,)
#Bb$class<-"intermediate"


#Cc<-as.data.frame(tidy(modelcont.phys,robust = TRUE))
#Cc<-Cc %>% "["(.,2:5,)
#Cc$class<-"physiological"
#Dd<-rbind(Aa,Bb,Cc)

#library(ggstance)
#pd=position_dodgev(height=0.3)
#ggplot(Dd,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="black")


#modelbin.funct <- brm(hyst.funct~ cent_pol+seed_cent+cent_minP+cent_flo_month +(1|name), data = dater, 
      #                 family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

#modelbin.phys <- brm(hyst.phys~ cent_pol+seed_cent+cent_minP+cent_flo_month+(1|name), data = dater, 
  #                    family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

#modelbin.inter <- brm(hyst.inter~ cent_pol+seed_cent+cent_minP+cent_flo_month+(1|name), data = dater, 
#                       family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000)


#aA<-as.data.frame(tidy(modelbin.funct,robust = TRUE))
#aA<-aA %>% "["(.,2:5,)
#aA$class<-"functional"

#bB<-as.data.frame(tidy(modelbin.inter,robust = TRUE))
#bB<-bB %>% "["(.,2:5,)
#bB$class<-"intermediate"


#cC<-as.data.frame(tidy(modelbin.phys,robust = TRUE))
#cC<-cC %>% "["(.,2:5,)
#cC$class<-"physiological"
#dD<-rbind(aA,bB,cC)

#ggplot(dD,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="black")




######interactions dont really need them
inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
colnames(dater)

modelcont.funct.int <- brm(offset.funct~ cent_pol+cent_floday+cent_minP+cent_floday:cent_pol+cent_floday:cent_minP+cent_pol:cent_minP+(1|name), data = dater, 
                     family = gaussian(), cov_ranef = list(name= A),iter=3000) 

#summary(modelcont.funct.int)


modelcont.phys.int <- brm(offset.phys~ cent_pol+cent_floday+cent_minP+cent_floday:cent_pol+cent_floday:cent_minP+cent_pol:cent_minP+(1|name), data = dater, 
                        family = gaussian(), cov_ranef = list(name= A),iter=3000) 

#summary(modelcont.phys.int)

#modelcont.inter.2<- brm(offset.inter~cent_pol+cent_floday+cent_minP+cent_floday:cent_pol+cent_floday:cent_minP+cent_pol:cent_minP+(1|name), data = dater, 
 #                        family = gaussian(), cov_ranef = list(name= A),iter=3000) 
#summary(modelcont.inter.2)

AAA<-as.data.frame(tidy(modelcont.funct.int,robust = TRUE))
AAA<-AAA %>% "["(.,1:7,)
AAA$class<-"functional"


CCC<-as.data.frame(tidy(modelcont.phys.int,robust = TRUE))
CCC<-CCC %>% "["(.,1:7,)
CCC$class<-"physiological"
DDD<-rbind(AAA,CCC)

DDD$term[which(DDD$term=="b_cent_pol")] <- "main effect: pollination syndrome"
DDD$term[which(DDD$term=="b_cent_floday")] <- "main effect: flowering time"
DDD$term[which(DDD$term=="b_cent_minP")] <- "main effect: minimum precipitation"
DDD$term[which(DDD$term=="b_cent_floday:cent_minP")] <- "interaction: flowering x precip."
DDD$term[which(DDD$term=="b_cent_pol:cent_minP")] <- "interaction: pollination x precip."
DDD$term[which(DDD$term=="b_cent_pol:cent_floday")] <- "interaction: pollination x flowering"


pd=position_dodgev(height=0.4)
HF.cont<-ggplot(DDD,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd,width=0.4)+geom_vline(aes(xintercept=0),color="black")+theme_bw()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("Continuous")

modelbin.funct.int <- brm(hyst.funct~ cent_pol+cent_floday+cent_minP+cent_floday:cent_pol+cent_floday:cent_minP+cent_pol:cent_minP+(1|name), data = dater, 
                         family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

modelbin.phys.int <- brm(hyst.phys~ cent_pol+cent_floday+cent_minP+cent_floday:cent_pol+cent_floday:cent_minP+cent_pol:cent_minP+(1|name), data = dater, 
                        family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

#modelbin.inter.2<- brm(hyst.inter~cent_pol+cent_floday+cent_minP+cent_floday:cent_pol+cent_floday:cent_minP+cent_pol:cent_minP+(1|name), data = dater, 
 #                       family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

aaa<-as.data.frame(tidy(modelbin.funct.int,robust = TRUE))
aaa<-aaa %>% "["(.,1:7,)
aaa$class<-"functional"

#bbb<-as.data.frame(tidy(modelbin.inter.2,robust = TRUE))
#bbb<-bbb %>% "["(.,1:7,)
#bbb$class<-"intermediate"


ccc<-as.data.frame(tidy(modelbin.phys.int,robust = TRUE))
ccc<-ccc %>% "["(.,1:7,)
ccc$class<-"physiological"
ddd<-rbind(aaa,ccc)

ddd$term[which(ddd$term=="b_cent_pol")] <- "main effect: pollination syndrome"
ddd$term[which(ddd$term=="b_cent_floday")] <- "main effect: flowering time"
ddd$term[which(ddd$term=="b_cent_minP")] <- "main effect: minimum precipitation"
ddd$term[which(ddd$term=="b_cent_floday:cent_minP")] <- "interaction: flowering x precip."
ddd$term[which(ddd$term=="b_cent_pol:cent_minP")] <- "interaction: pollination x precip."
ddd$term[which(ddd$term=="b_cent_pol:cent_floday")] <- "interaction: pollination x flowering"

HF.bin<-ggplot(ddd,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="black")+theme_bw()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("Binary")
library(gridExtra)
grid.arrange(HF.cont,HF.bin)

###############Sans interaction################################################################
modelcont.funct.noint <- brm(offset.funct~ cent_pol+cent_floday+cent_minP+(1|name), data = dater, 
                         family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelcont.phys.noint <- brm(offset.phys~ cent_pol+cent_floday+cent_minP+(1|name), data = dater, 
                        family = gaussian(), cov_ranef = list(name= A),iter=3000) 


#modelcont.inter.3<- brm(offset.inter~cent_pol+cent_floday+cent_minP+(1|name), data = dater, 
 #                       family = gaussian(), cov_ranef = list(name= A),iter=3000) 

summary(modelcont.funct.noint)
summary(modelcont.phys.noint)

extract_coefs<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,.25,.75,0.975))),"trait")
}


modelbin.funct.noint <- brm(hyst.funct~ cent_pol+cent_floday+cent_minP+(1|name), data = dater, 
                        family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

#summary(modelbin.funct.3)
modelbin.phys.noint <- brm(hyst.phys~ cent_pol+cent_floday+cent_minP+(1|name), data = dater, 
                      family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 


#modelbin.inter.3<- brm(hyst.inter~cent_pol+cent_floday+cent_minP+(1|name), data = dater, 
 #                       family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 


funct.cont<-extract_coefs(modelcont.funct.noint)
funct.cont$class<-"functional"
funct.bin<-extract_coefs(modelbin.funct.noint)
funct.bin$class<-"functional"

phys.cont<-extract_coefs(modelcont.phys.noint)
phys.bin<-extract_coefs(modelbin.phys.noint)
phys.cont$class<-"physiological"
phys.bin$class<-"physiological"


cont.noint<-rbind(phys.cont,funct.cont)
bin.noint<-rbind(phys.bin,funct.bin)



pd=position_dodgev(height=0.4)
hf.cont.noint<-ggplot(cont.noint,aes(Estimate,trait))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=class),linetype="solid",position=pd,width=0)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=class),linetype="dotted",position=pd,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("HF-Continuous")
hf.bin.noint<-ggplot(bin.noint,aes(Estimate,trait))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=class),linetype="solid",position=pd,width=0)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=class),linetype="dotted",position=pd,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("HF-Binary")


grid.arrange(hf.cont.noint,hf.bin.noint,nrow=1)

pp_check(modelcont.phys) ###models seem good
pp_check(modelcont.inter)
pp_check(modelbin.phys)
save.image("continuous.mods.RData")



stop("not an error, jsut excess below this")######
plot(dater$pol,dater$ave_flo_month)

AAAA<-as.data.frame(tidy(modelcont.funct.noint,robust = TRUE))
AAAA<-AAAA %>% "["(.,1:4,)
AAAA$class<-"functional"

#BBBB<-as.data.frame(tidy(modelcont.inter.3,robust = TRUE))
#BBBB<-BBBB %>% "["(.,1:4,)
#BBBB$class<-"intermediate"

CCCC<-as.data.frame(tidy(modelcont.phys.noint,robust = TRUE))
CCCC<-CCCC %>% "["(.,1:4,)
CCCC$class<-"physiological"
DDDD<-rbind(AAAA,CCCC)

ggplot(DDDD,aes(estimate,term))+geom_point(aes(color=class),position=pd)+geom_errorbarh(aes(xmin=lower,xmax=upper,color=class),position=pd)+geom_vline(aes(xintercept=0),color="black")

AAAAA<-as.data.frame(tidy(modelbin.funct.noint,robust = TRUE))
AAAAA<-AAAAA %>% "["(.,1:4,)
AAAAA$class<-"functional"

#BBBBB<-as.data.frame(tidy(modelbin.inter.3,robust = TRUE))
#BBBBB<-BBBBB %>% "["(.,1:4,)
#BBBBB$class<-"intermediate"
CCCCC<-as.data.frame(tidy(modelbin.phys.noint,robust = TRUE))
CCCCC<-CCCCC %>% "["(.,1:4,)
CCCCC$class<-"physiological"
DDDDD<-rbind(AAAAA,CCCCC)

