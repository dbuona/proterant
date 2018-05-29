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

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
mich.data$dev.time<-NA
mich.data$dev.time<-mich.data$fruiting-mich.data$flo_time

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
?phylo.d()
print(PhyloPro2)
#d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloPro3<-phylo.d(d,binvar=pro3)
PhyloPro3

###phlosignal for continuous trait
?phylosig()
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

#Mich.cent.funct<-phyloglm(pro2~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#                          start.beta=NULL, start.alpha=NULL,
 #                         boot=599,full.matrix = TRUE)


#summary(Mich.cent.funct)
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

#### Plot this
#bootest<-as.data.frame(Mich.cent.funct$coefficients)
#bootconf<-as.data.frame(Mich.cent.funct$bootconfint95)
#bootconf<-as.data.frame(t(bootconf))

#bootest<-rownames_to_column(bootest, "trait")
#bootconf<-rownames_to_column(bootconf, "trait")
#bootmich<-full_join(bootconf,bootest, by="trait")
#colnames(bootmich)<-c("trait","low","high","estimate")
#bootmich<-dplyr::filter(bootmich, trait!="alpha")
#bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
#bootmich$trait[bootmich$trait=="shade_bin"]<-"shade tolerance"
#bootmich$trait[bootmich$trait=="pol"]<-"pollination syndrome"
#bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
#bootmich$trait[bootmich$trait=="fruit_cent"]<-"fruit timing"
#bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"
#functplot<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
#functplot
############ main analysis plot using seed development as a predictor
cent.funct.seed<-phyloglm(pro2~pol+height_cent+flo_cent+dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)
summary(cent.funct.seed)
bootest<-as.data.frame(cent.funct.seed$coefficients)
bootconf<-as.data.frame(cent.funct.seed$bootconfint95)
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
bootmich$trait[bootmich$trait=="dev_time_cent"]<-"seed development"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

#jpeg("funct.effect.fill.jpeg")
functplot1<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
functplot1
#dev.off()
#####model 2 intermediate

#Mich.cent.intermed<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#                          start.beta=NULL, start.alpha=NULL,
#                          boot=599,full.matrix = TRUE)

#summary(Mich.cent.intermed)
#### Now Plot this
#bootest1<-as.data.frame(Mich.cent.intermed$coefficients)
#bootconf1<-as.data.frame(Mich.cent.intermed$bootconfint95)
#bootconf1<-as.data.frame(t(bootconf1))

#bootest1<-rownames_to_column(bootest1, "trait")
#bootconf1<-rownames_to_column(bootconf1, "trait")
#bootmich1<-full_join(bootconf1,bootest1, by="trait")
#colnames(bootmich1)<-c("trait","low","high","estimate")
#bootmich1<-dplyr::filter(bootmich1, trait!="alpha")
#bootmich1<-dplyr::filter(bootmich1, trait!="(Intercept)")
###names
#bootmich1$trait[bootmich1$trait=="shade_bin"]<-"shade tolerance"
#bootmich1$trait[bootmich1$trait=="pol"]<-"pollination syndrome"
#bootmich1$trait[bootmich1$trait=="height_cent"]<-"max height"
#bootmich1$trait[bootmich1$trait=="fruit_cent"]<-"fruit timing"
#bootmich1$trait[bootmich1$trait=="flo_cent"]<-"flower timing"

#interplot<-ggplot(bootmich1,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+theme(axis.text = element_text(size=14, hjust = .5))+xlim(-7,5)+guides(color="none")+ggtitle("Intermediate catagory")
#interplot

cent.intermed.seed<-phyloglm(pro~pol+height_cent+flo_cent+dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                            start.beta=NULL, start.alpha=NULL,
                            boot=599,full.matrix = TRUE)
summary(cent.intermed.seed)
bootest1<-as.data.frame(cent.intermed.seed$coefficients)
bootconf1<-as.data.frame(cent.intermed.seed$bootconfint95)
bootconf1<-as.data.frame(t(bootconf1))

bootest1<-rownames_to_column(bootest1, "trait")
bootconf1<-rownames_to_column(bootconf1, "trait")
bootmich1<-full_join(bootconf1,bootest1, by="trait")
colnames(bootmich1)<-c("trait","low","high","estimate")
bootmich1<-dplyr::filter(bootmich1, trait!="alpha")
bootmich1<-dplyr::filter(bootmich1, trait!="(Intercept)")
###names
bootmich1$trait[bootmich1$trait=="shade_bin"]<-"shade tolerance"
bootmich1$trait[bootmich1$trait=="pol"]<-"pollination syndrome"
bootmich1$trait[bootmich1$trait=="height_cent"]<-"max height"
bootmich1$trait[bootmich1$trait=="dev_time_cent"]<-"seed development"
bootmich1$trait[bootmich1$trait=="flo_cent"]<-"flower timing"

interplot1<-ggplot(bootmich1,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+theme(axis.text = element_text(size=14, hjust = .5))+xlim(-7,5)+guides(color="none")
interplot1


######model 3
#pro3 is response

#Mich.cent.phys<-phyloglm(pro3~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#                        start.beta=NULL, start.alpha=NULL,
#                        boot=599,full.matrix = TRUE)
#summary(Mich.cent.phys)
###and now this:
#bootest2<-as.data.frame(Mich.cent.phys$coefficients)
#bootconf2<-as.data.frame(Mich.cent.phys$bootconfint95)
#bootconf2<-as.data.frame(t(bootconf2))#

#bootest2<-rownames_to_column(bootest2, "trait")
#bootconf2<-rownames_to_column(bootconf2, "trait")
#bootmich2<-full_join(bootconf2,bootest2, by="trait")
#colnames(bootmich2)<-c("trait","low","high","estimate")
#bootmich2<-dplyr::filter(bootmich2, trait!="alpha")
#bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")
####names
#bootmich2$trait[bootmich2$trait=="shade_bin"]<-"shade tolerance"
#bootmich2$trait[bootmich2$trait=="pol"]<-"pollination syndrome"
#bootmich2$trait[bootmich2$trait=="height_cent"]<-"max height"
#bootmich2$trait[bootmich2$trait=="fruit_cent"]<-"fruit timing"
#bootmich2$trait[bootmich2$trait=="flo_cent"]<-"flower timing"

#physplot<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
#physplot

cent.phys.seed<-phyloglm(pro3~pol+height_cent+flo_cent+dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)
summary(cent.phys.seed)

bootest2<-as.data.frame(cent.phys.seed$coefficients)
bootconf2<-as.data.frame(cent.phys.seed$bootconfint95)
bootconf2<-as.data.frame(t(bootconf2))

bootest2<-rownames_to_column(bootest2, "trait")
bootconf2<-rownames_to_column(bootconf2, "trait")
bootmich2<-full_join(bootconf2,bootest2, by="trait")
colnames(bootmich2)<-c("trait","low","high","estimate")
bootmich2<-dplyr::filter(bootmich2, trait!="alpha")
bootmich2<-dplyr::filter(bootmich2, trait!="(Intercept)")
###names
bootmich2$trait[bootmich2$trait=="shade_bin"]<-"shade tolerance"
bootmich2$trait[bootmich2$trait=="pol"]<-"pollination syndrome"
bootmich2$trait[bootmich2$trait=="height_cent"]<-"max height"
bootmich2$trait[bootmich2$trait=="dev_time_cent"]<-"seed development"
bootmich2$trait[bootmich2$trait=="flo_cent"]<-"flower timing"
library(ggplot2)
jpeg("phys_effect_fill.jpeg")
physplot1<-ggplot(bootmich2,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+xlim(-7,5)+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
physplot1
dev.off()
plotty.all<-gridExtra::grid.arrange(functplot1, interplot1,physplot1,functplot,interplot,physplot, ncol=3,nrow=2)

############

##Unscale paramenters. You can rescale them by dividing betas by 2*sd corresponding x
    ##more work, but might allow for more reasonable average predictive comps
Mich.funct<-phyloglm(pro2~pol+heigh_height+flo_time+dev.time+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)
M2<-arm::standardize(Mich.funct)

summary(M2)
summary(Mich.funct)
###this is one way to rescale estimates
scaledestimates<-c(coef(Mich.funct)[2]/2*sd(mich.data$pol),
coef(Mich.funct)[3]/2*sd(mich.data$heigh_height),
coef(Mich.funct)[4]/2*sd(mich.data$flo_time),
coef(Mich.funct)[5]/2*sd(mich.data$dev.time),
coef(Mich.funct)[6]/2*sd(mich.data$shade_bin))
scale<-as.data.frame(scaledestimates)

scaledconf<-c(Mich.funct$bootconfint95[3]/2*sd(mich.data$pol),
Mich.funct$bootconfint95[4]/2*sd(mich.data$pol),
Mich.funct$bootconfint95[5]/2*sd(mich.data$heigh_height),
Mich.funct$bootconfint95[6]/2*sd(mich.data$heigh_height),
Mich.funct$bootconfint95[7]/2*sd(mich.data$flo_time),
Mich.funct$bootconfint95[8]/2*sd(mich.data$flo_time),
Mich.funct$bootconfint95[9]/2*sd(mich.data$dev.time),
Mich.funct$bootconfint95[10]/2*sd(mich.data$dev.time),
Mich.funct$bootconfint95[11]/2*sd(mich.data$shade_bin),
Mich.funct$bootconfint95[12]/2*sd(mich.data$shade_bin))
scale2<-as.data.frame(scaledconf)
scale2$pred<-c("pol","pol","heigh_height","heigh_height","flo_time","flo_time","dev.time","dev.time","shade_bin","shade_bin")
scale2$cont<-c("low","high","low","high","low","high","low","high","low","high")
scale3<-spread(scale2,cont,scaledconf)
scale<-rownames_to_column(scale, "pred")
goober<-full_join(scale,scale3, by="pred")
ggplot(goober,aes(scaledestimates,pred))+geom_point(size=2.5)+geom_segment(aes(y=pred,yend=pred,x=low,xend=high))+xlim(-1,1)+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")

###I don't think that was totally worth it.

#######side bar#################################################################
#Does pollination syndrome predict early flowering?

#a<-phylolm(flo_time~pol+fruit_cent+pro2+height_cent+shade_bin,data=mich.data,phy=mich.tree)
#summary(a)
#b<-phylolm(flo_time~pol+fruit_cent+pro+height_cent+shade_bin,data=mich.data,phy=mich.tree)
#summary(b)
#c<-phylolm(flo_time~pol+fruit_cent+pro3+height_cent+shade_bin,data=mich.data,phy=mich.tree)
#summary(c)
#d<-phylolm(flo_time~pol*pro,data=mich.data,phy=mich.tree)
#summary(d)
###Alt approach:
#Mich.cent.funct_winter<-phyloglm(pro2~pol*flo_cent+height_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#                          start.beta=NULL, start.alpha=NULL,
#                          boot=599,full.matrix = TRUE)
#summary(Mich.cent.funct_winter)

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
