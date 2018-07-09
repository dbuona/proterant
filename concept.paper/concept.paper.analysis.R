####c### This is the final analysis file for hysteranthy anaylsis on MTSV as of 3/28/18.
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

#if you dont want to run the model: 
load("hystmodels.RData")
#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time
###one more cleaninging tax
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)

mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)

#pro<- hysteranthy= before, and before with leaves 
  #phyloD =0.18 with 49 species hysteranthous

#pro2<-hysteranthy= before, before/with and with
    #phylod=0.06 with 101 species hysteranthous
#pro3<- hysteranthy=before only
    #D= 0.34 34 sp hysteranthous
###phylo.D###### calculate phylo d.
set.seed(122)
mich.tree$node.label<-NULL
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloD <- phylo.d(d, binvar=pro) ###regular hysteranthy
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

###Rescale predictors
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))
mich.data$dev_time_cent<-(mich.data$dev.time-mean(mich.data$dev.time))/(2*sd(mich.data$dev.time))
mich.data$tol_cent<-(mich.data$shade_bin-mean(mich.data$shade_bin))/(2*sd(mich.data$shade_bin))

###prepare for modeling
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
mich.data$height10<-mich.data$heigh_height/10
###models:


#####Mich.cent.funct uses pro2
##cent.funct.seed uses seed maturation duration rather than dispersal time
######Mich.cent.interm uses pro
####Mich.cent.phys
###Mih.fun

#number of bootstraps from Wilcox, R. R. (2010). Fundamentals of modern statistical methods: Substantially improving power and accuracy. Springer.


######do fruit and flo time covary? If yes removing one from the model should greatly change the effect
#test1<-phyloglm(pro2~pol+height_cent+flo_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
 #                         start.beta=NULL, start.alpha=NULL,
  #                        boot=599,full.matrix = TRUE)
#summary(test1)

#test2<-phyloglm(pro2~pol+height_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
 #               start.beta=NULL, start.alpha=NULL,
  #              boot=599,full.matrix = TRUE)

#summary(test2)
####The effect sizes remain relatively stable there fore we assume there is not a collinearity issue

############ main analysis plot using seed development as a predictor

z.funct.seed<-phyloglm(pro2~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

summary(z.funct.seed)

bootest2<-as.data.frame(z.funct.seed$coefficients)
bootconf2<-as.data.frame(z.funct.seed$bootconfint95)
bootconf2<-as.data.frame(t(bootconf2))

bootest2<-rownames_to_column(bootest2, "trait")
bootconf2<-rownames_to_column(bootconf2, "trait")
bootmich2<-full_join(bootconf2,bootest2, by="trait")
colnames(bootmich2)<-c("trait","low","high","estimate")
bootmich2<-dplyr::filter(bootmich2, trait!="alpha")
bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")
###names
bootmich2$trait[bootmich2$trait=="tol_cent"]<-"shade tolerance"
bootmich2$trait[bootmich2$trait=="pol_cent"]<-"pollination syndrome"
bootmich2$trait[bootmich2$trait=="height_cent"]<-"max height"
bootmich2$trait[bootmich2$trait=="dev_time_cent"]<-"seed development"
bootmich2$trait[bootmich2$trait=="flo_cent"]<-"flower timing"
bootmich2$class<-"functional-MTSV"

functplot1<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
functplot1

####
z.funct.seed.winter<-phyloglm(pro2~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent+pol_cent:flo_cent+dev_time_cent:flo_cent+height_cent:flo_cent+tol_cent:flo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)
summary(z.funct.seed.winter)
bootest<-as.data.frame(z.funct.seed.winter$coefficients)
bootconf<-as.data.frame(z.funct.seed.winter$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="alpha")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
bootmich$trait[bootmich$trait=="tol_cent"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol_cent"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="dev_time_cent"]<-"seed development"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"
bootmich$trait[bootmich$trait=="pol_cent:flo_cent"]<-"flower phenology x syndrome "
bootmich$trait[bootmich$trait=="flo_cent:dev_time_cent"]<-"flower phenology x development time"
bootmich$trait[bootmich$trait=="height_cent:flo_cent"]<-"flower phenology x height"
bootmich$trait[bootmich$trait=="flo_cent:tol_cent"]<-"flower phenology x tolerance"
bootmich$class<-"functional-MTSV"
functplot<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=1.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,7)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")+ggthemes::theme_few()
functplot

####marginal effects of flowering time on pollination

vels<-seq(0,1,1)
slopes <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[7]*vels
cbind(vels,slopes)
plot(vels, slopes, type = "l", lty = 1, ylim = c(-5, 0), xlab = "syndrome", ylab = "Marginal Effect of flower time")

###height
z.funct.seed.winter$coefficients
range(mich.data$height_cent)
height<-seq(-0.5603846, 1.0484774,0.1)
slopes2 <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[9]*height
cbind(height,slopes2)
plot(height, slopes2, type = "l", lty = 1, ylim = c(-7, 0), xlab = "height centered", ylab = "Marginal Effect of flower time")


###shade
z.funct.seed.winter$coefficients

shade<-seq(0,1,1)
slopes3 <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[10]*shade
cbind(shade,slopes3)
plot(shade, slopes3, type = "l", lty = 1, ylim = c(-4, 0), xlab = "shade centered", ylab = "Marginal Effect of flower time")

###seed development
z.funct.seed.winter$coefficients
range(mich.data$dev_time_cent)
seed<-seq(-0.7,2.1,0.2)
slopes4 <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[8]*seed
cbind(seed,slopes4)
plot(seed, slopes4, type = "l", lty = 1, ylim = c(-7, 2), xlab = "seed centered", ylab = "Marginal Effect of flower time")


###More plots for understanding interactions
interaction.plot(mich.data$pol,mich.data$flo_time,mich.data$pro2)
interaction.plot(mich.data$shade_bin,mich.data$flo_time,mich.data$pro2)

interaction.plot(mich.data$class,mich.data$flo_time,mich.data$pro2)

qplot(x = pol, y = pro2, data = mich.data, color = as.factor(flo_time)) +
  geom_smooth(method = "lm",se=FALSE), se=FALSE)+ggtitle("flower time X syndrome") 

qplot(x = shade_bin, y = pro2, data = mich.data, color = as.factor(flo_time)) +
  geom_smooth(method = "lm",se=FALSE)+ggtitle("flower time X shade tolerance") 

qplot(x = dev.time, y = pro2, data = mich.data, color = as.factor(flo_time)) +
  geom_smooth(method = "lm",se=FALSE)+ggtitle("flower time X development time") 

qplot(x = heigh_height ,y = pro2, data = mich.data, color = as.factor(flo_time)) +geom_smooth(method = "lm",se=FALSE)+ggtitle("flower time X height") 

###residuals and stuff

resid1<-as.data.frame(z.funct.seed.winter$fitted.values)
resid1<-rownames_to_column(resid1, "name")

resid2<-as.data.frame(z.funct.seed.winter$residuals)

resid2<-rownames_to_column(resid2, "name")
resid2$residuals<-resid$z.funct.seed.winter$residual
resid<-left_join(resid1,resid2)
qplot(z.funct.seed.winter$fitted.values,z.funct.seed.winter$residuals,data=resid)+geom_point()

goober<-filter(resid2,z.funct.seed.winter$residuals>=0.7)
goober2<-filter(resid2,z.funct.seed.winter$residuals<=-0.7)

residual.list<-rbind(goober,goober2)
write.csv(residual.list,"extreme.resid.csv")
##############################################################
z.phys.seed<-phyloglm(pro3~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)
summary(z.phys.seed)

bootest3<-as.data.frame(z.phys.seed$coefficients)
bootconf3<-as.data.frame(z.phys.seed$bootconfint95)
bootconf3<-as.data.frame(t(bootconf3))


bootest3<-rownames_to_column(bootest3, "trait")
bootconf3<-rownames_to_column(bootconf3, "trait")
bootmich3<-full_join(bootconf3,bootest3, by="trait")
colnames(bootmich3)<-c("trait","low","high","estimate")
bootmich3<-dplyr::filter(bootmich3, trait!="alpha")
bootmich3<-dplyr::filter(bootmich3, trait!="(Intercept)")
###names
bootmich3$trait[bootmich3$trait=="tol_cent"]<-"shade tolerance"
bootmich3$trait[bootmich3$trait=="pol_cent"]<-"pollination syndrome"
bootmich3$trait[bootmich3$trait=="height_cent"]<-"max height"
bootmich3$trait[bootmich3$trait=="dev_time_cent"]<-"seed development"
bootmich3$trait[bootmich3$trait=="flo_cent"]<-"flower timing"
bootmich3$class<-"physiological-MTSV"


z.phys.seed.winter<-phyloglm(pro3~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent+pol_cent:flo_cent+dev_time_cent:flo_cent+height_cent:flo_cent+tol_cent:flo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                           start.beta=NULL, start.alpha=NULL,
                           boot=599,full.matrix = TRUE)
summary(z.phys.seed.winter)
bootest1<-as.data.frame(z.phys.seed.winter$coefficients)
bootconf1<-as.data.frame(z.phys.seed.winter$bootconfint95)
bootconf1<-as.data.frame(t(bootconf1))


bootest1<-rownames_to_column(bootest1, "trait")
bootconf1<-rownames_to_column(bootconf1, "trait")
bootmich1<-full_join(bootconf1,bootest1, by="trait")
colnames(bootmich1)<-c("trait","low","high","estimate")
bootmich1<-dplyr::filter(bootmich1, trait!="alpha")
bootmich1<-dplyr::filter(bootmich1, trait!="(Intercept)")
###names
bootmich1$trait[bootmich1$trait=="tol_cent"]<-"shade tolerance"
bootmich1$trait[bootmich1$trait=="pol_cent"]<-"pollination syndrome"
bootmich1$trait[bootmich1$trait=="height_cent"]<-"max height"
bootmich1$trait[bootmich1$trait=="dev_time_cent"]<-"seed development"
bootmich1$trait[bootmich1$trait=="flo_cent"]<-"flower timing"
bootmich1$trait[bootmich1$trait=="pol_cent:flo_cent"]<-"flower phenology x syndrome "
bootmich1$trait[bootmich1$trait=="flo_cent:dev_time_cent"]<-"flower phenology x development time"
bootmich1$trait[bootmich1$trait=="height_cent:flo_cent"]<-"flower phenology x height"
bootmich1$trait[bootmich1$trait=="flo_cent:tol_cent"]<-"flower phenology x tolerance"
bootmich1$class<-"physiological-MTSV"

library(ggplot2)
#jpeg("phys_effect_fill.jpeg")
physplot1<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")+ggtitle("physical hyst")
physplot1
#dev.off()
###join the two for plotting togeteher
bootdata<-rbind(bootmich,bootmich1)
bootdata.nointer<-rbind(bootmich3,bootmich2)

library(ggstance)
pd=position_dodgev(height=0.25)
ggplot(bootdata,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()
ggplot(bootdata.nointer,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()

###add silvics back
silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")

####Silvics cleaning
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

###silve centering
silv.data$height_cent<-(silv.data$height-mean(silv.data$height))/(2*sd(silv.data$height))
silv.data$fruit_cent<-(silv.data$fruiting-mean(silv.data$fruiting))/(2*sd(silv.data$fruiting))
silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))
silv.data$pol_cent<-(silv.data$pol-mean(silv.data$pol))/(2*sd(silv.data$pol))
silv.data$tol_cent<-(silv.data$shade_bin-mean(silv.data$shade_bin))/(2*sd(silv.data$shade_bin))
silv.data$dev.time<-silv.data$fruiting-silv.data$flower_time
silv.data$dev_time_cent<-(silv.data$dev.time-mean(silv.data$dev.time))/(2*sd(silv.data$dev.time))
       
#####phylogenetic signals
e<-comparative.data(silv.tree,silv.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloSilv2<-phylo.d(e,binvar=pro2)
PhyloSilv2

PhyloSilv3<-phylo.d(e,binvar=pro3)
PhyloSilv3
silv.data<-  silv.data %>% remove_rownames %>% column_to_rownames(var="name")


###models
z.funct.silv<-phyloglm(pro2~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

summary(z.funct.silv)
bootestS<-as.data.frame(z.funct.silv$coefficients)
bootconfS<-as.data.frame(z.funct.silv$bootconfint95)
bootconfS<-as.data.frame(t(bootconfS))


bootestS<-rownames_to_column(bootestS, "trait")
bootconfS<-rownames_to_column(bootconfS, "trait")
bootsilvF<-full_join(bootconfS,bootestS, by="trait")
colnames(bootsilvF)<-c("trait","low","high","estimate")
bootsilvF<-dplyr::filter(bootsilvF, trait!="alpha")
bootsilvF<-dplyr::filter(bootsilvF, trait!="(Intercept)")
###names
bootsilvF$trait[bootsilvF$trait=="tol_cent"]<-"shade tolerance"
bootsilvF$trait[bootsilvF$trait=="pol_cent"]<-"pollination syndrome"
bootsilvF$trait[bootsilvF$trait=="height_cent"]<-"max height"
bootsilvF$trait[bootsilvF$trait=="dev_time_cent"]<-"seed development"
bootsilvF$trait[bootsilvF$trait=="flo_cent"]<-"flower timing"
bootsilvF$class<-"functional-Silvics"

z.phys.silv<-phyloglm(pro3~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,silv.data, silv.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

summary(z.phys.silv)
bootestS1<-as.data.frame(z.phys.silv$coefficients)
bootconfS1<-as.data.frame(z.phys.silv$bootconfint95)
bootconfS1<-as.data.frame(t(bootconfS1))


bootestS1<-rownames_to_column(bootestS1, "trait")
bootconfS1<-rownames_to_column(bootconfS1, "trait")
bootsilvP<-full_join(bootconfS1,bootestS1, by="trait")
colnames(bootsilvP)<-c("trait","low","high","estimate")
bootsilvP<-dplyr::filter(bootsilvP, trait!="alpha")
bootsilvP<-dplyr::filter(bootsilvP, trait!="(Intercept)")
###names
bootsilvP$trait[bootsilvP$trait=="tol_cent"]<-"shade tolerance"
bootsilvP$trait[bootsilvP$trait=="pol_cent"]<-"pollination syndrome"
bootsilvP$trait[bootsilvP$trait=="height_cent"]<-"max height"
bootsilvP$trait[bootsilvP$trait=="dev_time_cent"]<-"seed development"
bootsilvP$trait[bootsilvP$trait=="flo_cent"]<-"flower timing"
bootsilvP$class<-"physiological-Silvics"

bootdatasilv<-rbind(bootsilvF,bootsilvP)
pd=position_dodgev(height=0.25)

ggplot(bootdatasilv,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()

fullcomp<-rbind(bootdata.nointer,bootdatasilv)
pd=position_dodgev(height=0.4)
ggplot(fullcomp,aes(estimate,trait))+geom_point(size=2.5,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_bw()
###phyosig summary for 

(PhyloPro2) #0.06
PhyloPro3 #0.29
PhyloSilv2 #0.64
PhyloSilv3 #0.10
save.image(file="hystmodels.RData")

tab<-matrix(c(0.06,0.64,0.29,0.10))
rownames(tab)<-c("functional-MTSV","functional-Silvics","physiological-MTSV","Physiological-Silvics")
tab<-as.data.frame(tab)
tab<-rownames_to_column(tab, "Model")
###make this into a table at some point
###Which species are super early################################
superearl<-filter(mich.data, flo_time=="3.5")
superearl2<-filter(mich.data, flo_time=="4")
superearl3<-filter(mich.data, flo_time=="4.5")

stop("not an error, just ending the sourcing")


###########################################################################################

###average predictive comparsons(I dont think these are working yet.) Can I use divide by 4 rule?
coef(cent.funct.seed)
beta<-coef(cent.funct.seed)
hi<-1
lo<-0
###for [pollination syndrome]
beta[2]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$height_cent+beta[4]*mich.data$dev_time_cent+beta[5]*mich.data$fruit_cent+beta[6]*mich.data$shade_bin)-
    invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$height_cent+beta[4]*mich.data$flo_cent+beta[5]*mich.data$dev_time_cent+beta[6]*mich.data$shade_bin)

print(mean(delta))  #0.41


###for flowering
april<-(4-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
may<-(5-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
earl<-2*sd(mich.data$flo_time)*april+mean(mich.data$flo_time)
mid<-2*sd(mich.data$flo_time)*may+mean(mich.data$flo_time)   

delta2<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height_cent+beta[4]*april+beta[5]*mich.data$dev_time_cent+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height_cent+beta[4]*may+beta[5]*mich.data$dev_time_cent+beta[6]*mich.data$shade_bin)

mean(delta2) ##.24 which is right! yay

2*sd(mich.data$flo_time)*(-0.697978)+mean(mich.data$flo_time)
2*sd(mich.data$flo_time)*(-0.1677311)+mean(mich.data$flo_time)                                            
 
###how do I unscale this to make it meaningful this?
#=#z score = X-m/sd
#rescale SD(z) + M
#I did mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
#so I should do 2sd(z)+M


###APC unscaled#########################################

summary(Mich.funct)
coef(Mich.funct)
beta<-coef(Mich.funct)
hi<-1
lo<-0
###for [pollination syndrome]
beta[2]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$heigh_height+beta[4]*mich.data$flo_time+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$heigh_height+beta[4]*mich.data$flo_time+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)

mean(delta)  #0.40

###for flowering 
earl<-4
mid<-5
delta2<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$heigh_height+beta[4]*earl+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$heigh_height+beta[4]*mid+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)

mean(delta2)

#Whixh species are insect pollinated but not hysteranthous
buggy<-filter(mich.data, pol==0 & pro3==1)
unique(buggy$Family)
###tropical Rhamnaceae, Lauraceae, Anacardiaceae, Rutaceae
phy<-as.data.frame(coef(cent.funct.seed.winter))
phy<- phy %>%  rownames_to_column(var="effect")

nophy<-as.data.frame(coef(cent.funct.seed.nophylo))
nophy<- nophy %>%  rownames_to_column(var="effect")
phycom<-left_join(phy,nophy)

#-------------------------------------------------------#
#okay no to see if I can replicate the results of the scoop
###different hysteranthy definition
####different predictors
###didn't z-score (they log transformed continuous variables)
###reduced species list
#no interaactions
### they conflate early flowering hypotheses with hysteranthous hypotheses

### also I think they maybe ran a seperate model for each hypothesis

#1 ###take out flower time and compare functional to physiological
#predictor seed dev should be big and pol samall
funct.check<-phyloglm(pro2~pol+heigh_height+dev.time+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=10,full.matrix = TRUE)
summary(funct.check)
phys.check<-phyloglm(pro3~pol+heigh_height+dev.time+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=10,full.matrix = TRUE)
summary(phys.check)
##its true that dev time is increasing in importance and syndrome decreasing with different definitions
#### try it z scored
z.check.funct<-phyloglm(pro2~pol+height_cent++dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=10,full.matrix = TRUE)
summary(z.check.funct)
#pol 2.88, dev -0.5
z.check.phys<-phyloglm(pro3~pol+height_cent++dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=10,full.matrix = TRUE)
summary(z.check.phys)
#pol 2.01 dev -1.7

####lets look at dev time flower time and polinaton
z.check.funct2<-phyloglm(pro2~pol+flo_cent+dev_time_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=10,full.matrix = TRUE)

z.check.phys2<-phyloglm(pro3~pol+flo_cent+dev_time_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=10,full.matrix = TRUE)
summary(z.check.funct2)
##pol 2.69, flo, -2.98, dev time -0.9
summary(z.check.phys2) 
##pol 1.5, flo -4.11, dev -1.43

funct.check2<-phyloglm(pro2~pol+flo_time+dev.time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=10,full.matrix = TRUE)
summary(funct.check2)
#pol 2.7, flo -1.86 dev -0.18

phys.check2<-phyloglm(pro3~pol+flo_time+dev.time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=10,full.matrix = TRUE)
summary(phys.check2)
#pol 0.99, flo -2.46 dev -.22

####summary accross the board, polination syndrome decrease with phys definition

