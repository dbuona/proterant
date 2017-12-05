###This is a tidy analysis compiling and comparing all of Dan B's hysteranthy codes from elsewhere
###Began 5 Spet 2017

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

silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")
##Mich trees phylo D
set.seed(122)
mich.tree$node.label<-NULL
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloD <- phylo.d(d, binvar=pro) ###Physiological hysteranthy
PhyloD
##functionalhysteranthy
PhyloPro2<-phylo.d(d,binvar=pro2)
PhyloPro2
###Silvics phyloD
e<-comparative.data(silv.tree,silv.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloE <- phylo.d(e, binvar=pro)
PhyloE

#######Michigan cleaning######################
#clean av fruit time
mich.data$av_fruit_time[mich.data$av_fruit_time=="persistant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="persitant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="unreported"]<-9    
mich.data$av_fruit_time<-as.numeric(mich.data$av_fruit_time)

###clean fruiting
#mich.data$fruiting[mich.data$fruiting==19]<-7

#####add a new column for a adjusting for red acorn time
mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9                                      
mich.data$fruiting<-as.numeric(mich.data$fruiting)

####Silvics cleaning
###fruiting
silv.data$fruiting<-NA
silv.data$fruiting<-silv.data$av_fruit_time
silv.data$fruiting[silv.data$fruiting==21]<-9

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

###super biological hysreanthy
mich.data["pro3"]<-NA
mich.data$pro3[mich.data$Phen.sequence == "pro"] <- 1
mich.data$pro3[mich.data$Phen.sequence == "pro/syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn/ser"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "ser"] <- 0 
mich.data$pro3[mich.data$Phen.sequence== "hyst"] <- 0
##super biological for michigan trees
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloPro3<-phylo.d(d,binvar=pro3)
PhyloPro3


###functional phylo.D for suilvics
e<-comparative.data(silv.tree,silv.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloSilv2<-phylo.d(e,binvar=pro2)
PhyloSilv2

PhyloSilv3<-phylo.d(e,binvar=pro3)
PhyloSilv3

#####Centering
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))


silv.data$height_cent<-(silv.data$height-mean(silv.data$height))/(2*sd(silv.data$height))
silv.data$fruit_cent<-(silv.data$fruiting-mean(silv.data$fruiting))/(2*sd(silv.data$fruiting))
silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))




mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<-  silv.data %>% remove_rownames %>% column_to_rownames(var="name")


###uncentered
mich5<-phyloglm(pro~pol+heigh_height+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5)

#centering just height to  compare with silvics
mich.data$height10<-mich.data$heigh_height-mean(mich.data$heigh_height)
mich5h<-phyloglm(pro~pol+height10+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5h)

###average predictive comparisons (These are very sensative to adjusting)
beta<-coef(mich5h)
hi<-1
lo<-0
###for [pollination syndrome]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$height10+beta[4]*mich.data$flo_time+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$height10+beta[4]*mich.data$flo_time+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)
print(mean(delta))

###for flowering
earl<-4
mid<-5
delta<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height10+beta[4]*earl+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height10+beta[4]*mid+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)

print(mean(delta))


###centered for comparision between bianry and continuous
mich5cent<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(mich5cent)


###Functional hysteranthy-centered
Mich5cent.funct<-phyloglm(pro2~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)
summary(Mich5cent.funct)

Mich5cent.super<-phyloglm(pro3~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)
summary(Mich5cent.super)

###uncentered silvics doesn't run with normal height, subtract mean
silv.data$height10<-silv.data$height-mean(silv.data$height)
sil5h<-phyloglm(pro~pol+flower_time+height10+fruiting+shade_bin,silv.data, silv.tree, method = "logistic_MPLE", btol = 400, log.alpha.bound = 1,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(sil5h)
beta<-coef(sil5h)
hi<-1
lo<-0
###for [pollination syndrome]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*silv.data$flower_time+beta[4]*silv.data$height10+beta[5]*silv.data$fruiting+beta[6]*silv.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*silv.data$flower_time+beta[4]*silv.data$height10+beta[5]*silv.data$fruiting+beta[6]*silv.data$shade_bin)
print(mean(delta))

###for flowering
earl<-4
mid<-5

delta<-invlogit(beta[1]+beta[2]*silv.data$pol+beta[3]*earl+beta[4]*silv.data$height10+beta[5]*silv.data$fruiting+beta[6]*silv.data$shade_bin)-
  invlogit(beta[1]+beta[2]*silv.data$pol+beta[3]*mid+beta[4]*silv.data$height10+beta[5]*silv.data$fruiting+beta[6]*silv.data$shade_bin)
print(mean(delta))

#fruiting
summer<-8
winter<-11

delta<-invlogit(beta[1]+beta[2]*silv.data$pol+beta[3]*silv.data$flower_time+beta[4]*silv.data$height10+beta[5]*summer+beta[6]*silv.data$shade_bin)-
  invlogit(beta[1]+beta[2]*silv.data$pol+beta[3]*silv.data$flower_time+beta[4]*silv.data$height10+beta[5]*winter+beta[6]*silv.data$shade_bin)
print(mean(delta))


###centered silvics
sil5.cent<-phyloglm(pro~pol+flo_cent+height_cent+fruit_cent+shade_bin,silv.data, silv.tree, method = "logistic_MPLE", btol = 60, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot=50,full.matrix = TRUE)
summary(sil5.cent)
cor(silv.data$fruit_cent,silv.data$flo_cent)

#### functional silvics
sil5.cent.funct<-phyloglm(pro2~pol+flo_cent+height_cent+fruit_cent+shade_bin,silv.data, silv.tree, method = "logistic_MPLE", btol = 60, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(sil5.cent.funct)

sil5.cent.super<-phyloglm(pro3~pol+flo_cent+height_cent+fruit_cent+shade_bin,silv.data, silv.tree, method = "logistic_MPLE", btol = 60, log.alpha.bound = 4,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)
summary(sil5.cent.super)


###centered full model

##centered full model with oaks original
mich5centoaks<-phyloglm(pro~pol+height_cent+flo_cent+av_fruit_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(mich5centoaks)


#plotting.
###the tree and variable:
#plotting

library(gridExtra)
summary(mich5cent)

###Non-boot strapped
est<-as.data.frame(coef(mich5cent))
est<-rownames_to_column(est, "name")
ints<-as.data.frame(confint(mich5cent,level = 0.95))
ints<-rownames_to_column(ints, "name")
colnames(ints)[2] <- "low"
colnames(ints)[3] <- "high"
colnames(est)[2] <- "estimate"
foo<-left_join(est,ints)
foo<-filter(foo,estimate<10)
plotI<-ggplot(foo,aes(estimate,name))+geom_point()+geom_segment(aes(y=name,yend=name,x=low,xend=high))+ggtitle("Main effects of predictors on Hysteranthy: MTSV")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")

###plot for silvics
est2<-as.data.frame(coef(sil5.cent))
est2<-rownames_to_column(est2, "name")
ints2<-as.data.frame(confint(sil5.cent,level = 0.95))
ints2<-rownames_to_column(ints2, "name")
colnames(ints2)[2] <- "low"
colnames(ints2)[3] <- "high"
colnames(est2)[2] <- "estimate"
foo2<-left_join(est2,ints2)
foo2<-filter(foo2,estimate<10)
plotII<-ggplot(foo2,aes(estimate,name))+geom_point()+geom_segment(aes(y=name,yend=name,x=low,xend=high))+ggtitle("Main effects of predictors on Hysteranthy: Silvics")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")
grid.arrange(plotI,plotII, ncol=2)
###bootstrapped
bootest<-as.data.frame(mich5cent$coefficients)
bootconf<-as.data.frame(mich5cent$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "name")
bootconf<-rownames_to_column(bootconf, "name")
bootmich<-full_join(bootconf,bootest, by="name")
colnames(bootmich)<-c("name","low","high","estimate")
bootmich<-dplyr::filter(bootmich, name!="alpha")
bootmich<-dplyr::filter(bootmich, name!="(Intercept)")
###names
bootmich$name[bootmich$name=="shade_bin"]<-"shade tolerance"
bootmich$name[bootmich$name=="pol"]<-"pollination syndrome"
bootmich$name[bootmich$name=="height_cent"]<-"max height"
bootmich$name[bootmich$name=="fruit_cent"]<-"fruit timing"
bootmich$name[bootmich$name=="flo_cent"]<-"flower timing"

ggplot(bootmich,aes(estimate,name))+geom_point()+geom_segment(aes(y=name,yend=name,x=low,xend=high))+ggtitle("Main effects of predictors on Hysteranthy: MTSV")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")

#### do it for silvics
bootest<-as.data.frame(sil5.cent$coefficients)
bootconf<-as.data.frame(sil5.cent$bootconfint95)
bootconf<-as.data.frame(t(bootconf))
bootest<-rownames_to_column(bootest, "name")
bootconf<-rownames_to_column(bootconf, "name")
bootsil<-full_join(bootconf,bootest, by="name")
colnames(bootsil)<-c("name","low","high","estimate")
bootsil<-dplyr::filter(bootsil, name!="alpha")
bootsil<-dplyr::filter(bootsil, name!="(Intercept)")

bootsil$name[bootsil$name=="shade_bin"]<-"shade tolerance"
bootsil$name[bootsil$name=="pol"]<-"pollination syndrome"
bootsil$name[bootsil$name=="height_cent"]<-"max height"
bootsil$name[bootsil$name=="fruit_cent"]<-"fruit timing"
bootsil$name[bootsil$name=="flo_cent"]<-"flower timing"

ggplot(bootsil,aes(estimate,name))+geom_point()+geom_segment(aes(y=name,yend=name,x=low,xend=high))+ggtitle("Main effects of predictors on Hysteranthy: Silvics")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")



##This plotting is based on bootstrap cint but the data sheet is not updated
#boot<-read.csv("mich5bootoutput.csv",header=TRUE)
#head(boot)
#boot<-dplyr::filter(boot,Coefficients!="(Intercept)")
#ggplot(boot,aes(Estimate,Coefficients))+geom_point()+geom_segment(aes(y=Coefficients,yend=Coefficients,x=lowerbootCI,xend=upperbootCI))+ggtitle("Main effects of predictors on Hysteranthy with bootstrap")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")

####binomial plots

plot1<-ggplot(mich.data, aes(x=fruiting, y=pro)) + geom_point(shape=1, position=position_jitter(width=.05,height=.05))  + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

plot2<-ggplot(mich.data, aes(x=flo_time, y=pro)) + geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

plot3<-ggplot(mich.data, aes(x=pol, y=pro)) + 
  geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

plot4<-ggplot(mich.data, aes(x=shade_bin, y=pro)) + 
  geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

plot5<-ggplot(mich.data, aes(x=heigh_height, y=pro)) + 
  geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

grid.arrange(plot1, plot2,plot3,plot4,plot5, ncol=2)

###trying to plot hysteranthy as binary adn ugly af
z <- as.factor(mich.data$pro); names(z)<-mich.tree$tip.label
#create a vector to be filled with colors for plotting purposes, matching the diel activity pattern
mycol<-character(length(z))
mycol[mich.data$pro==1]<-"black"
mycol[mich.data$pro==0]<-"red"

#make very narrow margins
par(mar=c(0,0,0,0))

#and now plot the phylogeny and the info on diel activity pattern
plot(mich.tree,x.lim=c(0,400),show.tip.label=FALSE)
tiplabels(pch=22, col="black", bg=mycol, cex=.5)
