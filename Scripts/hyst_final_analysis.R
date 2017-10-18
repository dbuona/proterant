###This is a tidy analysis compiling and comparing all of Dan B's hysteranthy codes from elsewhere
###Began 5 Spet 2017

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(tidyverse)
library(boot)
library(phylolm)
library(ggplot2)

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")

silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")


#keeler.tree<-read.tree("pruned_keeler.tre")
#keeler.data<-read.csv("keeler_cleaned.csv")

#setdiff(keeler.data$name,mich.data$name)
###clean fruiting
mich.data$fruiting[mich.data$fruiting==19]<-7

#####add a new column for a adjusting for red acorn time
mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9                                      
mich.data$fruiting<-as.numeric(mich.data$fruiting)

silv.data$fruiting<-NA
silv.data$fruiting<-silv.data$av_fruit_time
silv.data$fruiting[silv.data$fruiting==21]<-9

#####Centering
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))

silv.data$height_cent<-(silv.data$height-mean(silv.data$height))/(2*sd(silv.data$height))
silv.data$fruit_cent<-(silv.data$fruiting-mean(silv.data$fruiting))/(2*sd(silv.data$fruiting))
silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))



###########compare 2 variable models####
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<-  silv.data %>% remove_rownames %>% column_to_rownames(var="name")

mich2<-phyloglm(pro~pol+flo_time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=100,full.matrix = TRUE)

sil2<-phyloglm(pro~pol+flower_time,silv.data, silv.tree, method = "logistic_MPLE", btol = 60, log.alpha.bound = 4,
               start.beta=NULL, start.alpha=NULL,
               boot=50,full.matrix = TRUE)
summary(sil2)

summary(mich2)
coef(sil2)
coef(mich2)
####pollination only become significant when av_fruit_time is in the model

mich5<-phyloglm(pro~pol+heigh_height+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5)
sil5<-phyloglm(pro~pol+flo_cent+height_cent+fruit_cent+shade_bin,silv.data, silv.tree, method = "logistic_MPLE", btol = 60, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot=50,full.matrix = TRUE)
summary(sil5)

coef(sil5)
coef(mich5cent)
coef(mich5)
###centered full model
mich5cent<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5cent)

###merge
colnames(mich.data)[which(names(mich.data) == "heigh_height")] <- "height"
colnames(mich.data)[which(names(mich.data) == "flo_time")] <- "flower_time"

mich.data$ds<-"michigan"
silv.data$ds<-"silvics"

mich.data<-dplyr::select(mich.data,pro,pol,height,height_cent,flower_time,fruiting,fruit_cent,flo_cent, shade_bin,ds)
silv.data<-dplyr::select(silv.data,pro,pol,height,height_cent,flower_time,fruiting,fruit_cent,flo_cent, shade_bin,ds)
mich.data<-rownames_to_column(mich.data, "name")
silv.data<-rownames_to_column(silv.data, "name")
bigdata<-rbind(mich.data,silv.data)

big1<-glm(pro~pol+height_cent+flower_time+fruiting+shade_bin+ds,family = binomial(link="logit"),data=bigdata)
summary(big1)

###can we do this phylogenetically?
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
names.intree<-treee$tip.label

# list of my species myspecies
namelist<-unique(bigdata$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
#plot(pruned.by.anthy)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #107 species include

addins<-setdiff(namelist,mytree.names) #30 species did not make it
###
###make ultrametric (using mean path length smoothing, could also try penalized maximum likelihood with chronos())
is.ultrametric(pruned.by.anthy)
help(chronoMPL)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)
#plot(pruned.by.anthy)
#adding species to tree at root
species<-addins
for(i in 1:length(species)) pruned.by.anthy<-add.species.to.genus(pruned.by.anthy,species[i],
                                                                  where="root")
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #107 species include
final.df<-bigdata[match(mytree.names, bigdata$name),]
namelist3<-final.df$name
namelist3==mytree.names
final.df$name== mytree.names

####it can seem to handele duplicate data
biggie<-phyloglm(pro~pol+height_cent+flower_time+fruiting+shade_bin+ds,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
         start.beta=NULL, start.alpha=NULL,
         boot=100,full.matrix = TRUE)


#plotting.
###the tree and variable:
#plotting

summary(mich5)

coef(mich5)
est<-as.data.frame(coef(mich5))
est<-rownames_to_column(est, "name")
ints<-as.data.frame(confint(mich5,level = 0.95))
ints<-rownames_to_column(ints, "name")
colnames(ints)[2] <- "low"
colnames(ints)[3] <- "high"
colnames(est)[2] <- "estimate"
foo<-left_join(est,ints)
foo<-filter(foo,estimate<10)
ggplot(foo,aes(estimate,name))+geom_point()+geom_segment(aes(y=name,yend=name,x=low,xend=high))+ggtitle("Main effects of predictors on Hysteranthy")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")

boot<-read.csv("mich5bootoutput.csv",header=TRUE)
head(boot)
boot<-dplyr::filter(boot,Coefficients!="(Intercept)")
ggplot(boot,aes(Estimate,Coefficients))+geom_point()+geom_segment(aes(y=Coefficients,yend=Coefficients,x=lowerbootCI,xend=upperbootCI))+ggtitle("Main effects of predictors on Hysteranthy with bootstrap")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+guides(color="none")

####binomial plots
ggplot(mich.data, aes(x=fruiting, y=pro)) + geom_point(shape=1, position=position_jitter(width=.05,height=.05))  + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

ggplot(mich.data, aes(x=flo_time, y=pro)) + geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

ggplot(mich.data, aes(x=pol, y=pro)) + 
  geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

ggplot(mich.data, aes(x=shade_bin, y=pro)) + 
  geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()

ggplot(mich.data, aes(x=heigh_height, y=pro)) + 
  geom_point(shape=1, position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)+theme_minimal()
