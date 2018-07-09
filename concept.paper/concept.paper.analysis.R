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
cent.funct.seed<-phyloglm(pro2~pol+height_cent+flo_cent+dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)
summary(cent.funct.seed)
cent.funct.seed.nophylo<-glm(pro2~pol+height_cent+flo_cent+dev_time_cent+shade_bin+flo_cent:pol, data=mich.data, family=binomial())
summary(cent.funct.seed.nophylo)
#### dev time is the only thing that changes (stronger in phylocorrected)
####s, I didn't use the standardized bianary predictor: THis is the model I should use, though nothing changes that much
z.funct.seed<-phyloglm(pro2~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

summary(z.funct.seed)
summary(cent.funct.seed)

bootest<-as.data.frame(z.funct.seed$coefficients)
bootconf<-as.data.frame(z.funct.seed$bootconfint95)
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

#jpeg("funct.effect.fill.jpeg")
functplot1<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
functplot1
#dev.off()
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

functplot<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=1.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,7)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")+ggthemes::theme_few()
functplot

#funct.seed.winter<-phyloglm(pro2~pol+heigh_height+flo_time+dev.time+shade_bin+pol:flo_time+dev.time:flo_time+heigh_height:flo_time+shade_bin:flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                              start.beta=NULL, start.alpha=NULL,
                              boot=59,full.matrix = TRUE)

summary(funct.seed.winter)
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

#cent.phys.seed<-phyloglm(pro3~pol+height_cent+flo_cent+dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)
summary(cent.phys.seed)

z.phys.seed.winter<-phyloglm(pro3~pol_cent+height_cent+flo_cent+dev_time_cent+tol_cent+pol_cent:flo_cent+dev_time_cent:flo_cent+height_cent:flo_cent+tol_cent:flo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                             start.beta=NULL, start.alpha=NULL,
                             boot=599,full.matrix = TRUE)
summary(z.phys.seed.winter)
bootest<-as.data.frame(z.phys.seed.winter$coefficients)
bootconf<-as.data.frame(z.phys.seed.winter$bootconfint95)
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


library(ggplot2)
#jpeg("phys_effect_fill.jpeg")
physplot1<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")+ggtitle("physical hyst")
physplot1
#dev.off()
plotty.all<-gridExtra::grid.arrange(functplot1, interplot1,physplot1,functplot,interplot,physplot, ncol=3,nrow=2)
############
work<-ls()
save(list =work, file="hystmodels.RData")

###Which species are super early
superearl<-filter(mich.data, flo_time=="3.5")
superearl2<-filter(mich.data, flo_time=="4")
superearl3<-filter(mich.data, flo_time=="4.5")

stop("not an error, just ending the sourcing")

####question should I bin flowering time different
mich.data$smooth_flo<-mich.data$flo_time

mich.data$smooth_flo[which(mich.data$flo_time==3.5 )] <- 3
mich.data$smooth_flo[which(mich.data$flo_time==4.5 )] <- 4
mich.data$smooth_flo[which(mich.data$flo_time==5.5 )] <- 5
mich.data$smooth_flo[which(mich.data$flo_time==6.5 )] <- 6
mich.data$smooth_flo[which(mich.data$flo_time==7.5 )] <- 7
mich.data$smooth_flo[which(mich.data$flo_time==10.5 )] <- 10

colnames(mich.data)
mich.data$smoothflo_cent<-(mich.data$smooth_flo-mean(mich.data$smooth_flo))/(2*sd(mich.data$smooth_flo))
z.funct.smooth.winter<-phyloglm(pro2~pol_cent+height_cent+smoothflo_cent+dev_time_cent+tol_cent+pol_cent:smoothflo_cent+height_cent:smoothflo_cent+dev_time_cent:smoothflo_cent+tol_cent:smoothflo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 40, log.alpha.bound = 4,
                              start.beta=NULL, start.alpha=NULL,
                              boot=599,full.matrix = TRUE)

summary(z.funct.smooth.winter)
bootest<-as.data.frame(z.funct.smooth.winter$coefficients)
bootconf<-as.data.frame(z.funct.smooth.winter$bootconfint95)
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
bootmich$trait[bootmich$trait=="smoothflo_cent"]<-"flower timing"
bootmich$trait[bootmich$trait=="pol_cent:smoothflo_cent"]<-"flower phenology x syndrome "
bootmich$trait[bootmich$trait=="smoothflo_cent:dev_time_cent"]<-"flower phenology x development time"
bootmich$trait[bootmich$trait=="height_cent:smoothflo_cent"]<-"flower phenology x height"
bootmich$trait[bootmich$trait=="smoothflo_cent:tol_cent"]<-"flower phenology x tolerance"
funcsmooth<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")+ggtitle("smooth hyst")
funcsmooth

####phsically smooth
z.inter.smooth.winter<-phyloglm(pro~pol_cent+height_cent+smoothflo_cent+dev_time_cent+tol_cent+pol_cent:smoothflo_cent+height_cent:smoothflo_cent+dev_time_cent:smoothflo_cent+tol_cent:smoothflo_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 40, log.alpha.bound = 4,
                                start.beta=NULL, start.alpha=NULL,
                                boot=599,full.matrix = TRUE)

summary(z.inter.smooth.winter)
bootest<-as.data.frame(z.inter.smooth.winter$coefficients)
bootconf<-as.data.frame(z.inter.smooth.winter$bootconfint95)
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
bootmich$trait[bootmich$trait=="smoothflo_cent"]<-"flower timing"
bootmich$trait[bootmich$trait=="pol_cent:smoothflo_cent"]<-"flower phenology x syndrome "
bootmich$trait[bootmich$trait=="smoothflo_cent:dev_time_cent"]<-"flower phenology x development time"
bootmich$trait[bootmich$trait=="height_cent:smoothflo_cent"]<-"flower phenology x height"
bootmich$trait[bootmich$trait=="smoothflo_cent:tol_cent"]<-"flower phenology x tolerance"
intersmooth<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")+ggtitle("smooth inter hyst")
intersmooth



save(list =work, file="hystmodels.RData")
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

