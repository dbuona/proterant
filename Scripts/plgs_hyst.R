###first attempt at plgs
## phylogeny of proterany as a bianary trait, adapt from first phylo_script.R april 10
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()
#load packages
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(phylobase)
library(tidyverse)
setwd("~/Documents/git/proterant/input")

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
#list of species in tree
names.intree<-treee$tip.label

### read in data from michigan trees, format it like Zanne
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy$name<-paste(anthy$Genus,anthy$Species,sep="_")
# list of my species myspecies
namelist<-unique(anthy$name)
##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

#put the tree and data sheet in same order
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names
##okay all catagories match
#now make a new column for proteranthy as a bianary trait
final.df["pro"]<-NA
final.df$pro[final.df$Phen.sequence == "pro"] <- 1
final.df$pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$pro[final.df$Phen.sequence== "syn"] <- 0
final.df$pro[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$pro[final.df$Phen.sequence== "ser"] <- 0 
final.df$pro[final.df$Phen.sequence== "hyst"] <- 0

###now make pollination syndrom bianary
##change acer
final.df$Pollination[final.df$name == "Acer_spicatum"] <- "insect"
final.df$Pollination[final.df$name == "Acer_pensylvanicum"] <- "insect"
final.df$Pollination[final.df$name == "Acer_negundo"] <- "wind"
final.df$Pollination[final.df$name == "Acer_platanoides"] <- "wind"
final.df$Pollination[final.df$name == "Acer_rubrum"] <- "wind"
final.df$Pollination[final.df$name == "Acer_saccharinum"] <- "wind"
final.df$Pollination[final.df$name == "Acer_saccharum"] <- "wind"

final.df["pol"]<-0
final.df$pol[final.df$Pollination == "insect"] <- 0
final.df$pol[final.df$Pollination == "wind"] <- 1

###now height
final.df$class<-NA
final.df<- within(final.df, class[heigh_height<10]<-1)
final.df<- within(final.df, class[heigh_height>=10 & heigh_height <20]<-2)
final.df<- within(final.df, class[heigh_height>=20]<-3)


###make factors for bianary traits
final.df$pro<-as.integer(final.df$pro)
final.df$pol<-as.integer(final.df$pol)
#final.df$class<-as.integer(final.df$class)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
head(final.df)
name.check(pruned.by.anthy, final.df)


library("phylolm")
##cleaning

##############################################################################################################
###models: for 1 predictor: pollination syndrone
##How to I interprest this model with catagorical data? DO I need a dummy varaible? possibly helkpful https://onlinecourses.science.psu.edu/stat504/node/225
#mod1<-glm(pro~Pollination-1,family = binomial(link="logit"),data=final.df)
#mod2<-phyloglm(pro~as.factor-1,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
              #start.beta=NULL, start.alpha=NULL,
              #boot = 0, full.matrix = TRUE)
#mod3<-phyloglm(pro~Pollination+class,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 20, log.alpha.bound = 4,
               #start.beta=NULL, start.alpha=NULL,
#####################################################################################
######model for heigh height only
mod4<-glm(pro~heigh_height,family = binomial(link="logit"),data=final.df)
summary(mod4)
inv.logit(coef(mod4)[1]+coef(mod4)[2]*mean(final.df$heigh_height))
mean(final.df$heigh_height)
###results:
#=.3420594 at mean height of 20.275 meters
#probability of switch to proteranthy with every increased meter in heigh_height
.06413/4
#0.016325

mod4a<-phyloglm(pro~heigh_height,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
start.beta=NULL, start.alpha=NULL,
boot = 0, full.matrix = TRUE)
summary(mod4a)
0.034187/4
#0.00854675

###now as catogorical
mod5<-glm(pro~class,family = binomial(link="logit"),data=final.df)
summary(mod5)
0.3788/4
#=0.09

mod5a<-phyloglm(pro~class,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod5a)
-0.00034306/4
#-8.5765e-05
##############bianaryxbianary for polliantion syndrome and hysteranthy##http://data.princeton.edu/wws509/notes/c3.pdf
mod6<-glm(pro~pol,family = binomial(link="logit"),data=final.df)
summary(mod6)
1.5311/4
##.382775 ~40%
mod6a<-phyloglm(pro~pol,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod6a)
0.78879/4
###0.197 ~20% chance you become hysteranthous if you switch from insect to wind pollination

##########################################################height and pollination syndrome#############
mod7<-glm(pro~pol+heigh_height,family = binomial(link="logit"),data=final.df)
summary(mod7)
1.35708/4 #=0.33927
0.02552/4 #=0.00638

mod7a<-phyloglm(pro~pol+heigh_height,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 30, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod7a)
1.049309/4  #=0.2623273
0.015686/4 #=0.0039215

#rescale heigh_height
height10<-final.df$heigh_height/10
mod8<-glm(pro~pol+height10,family = binomial(link="logit"),data=final.df)
summary(mod8)
1.3571/4 #=0.339275
0.2552/4 #0.0638
mod8a<-phyloglm(pro~pol+height10,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 30, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod8a)
0.47130/4 #=0.117825
0.31157/4 # 0.0778925

######interactions#####3
mod9<-glm(pro~pol*heigh_height,family = binomial(link="logit"),data=final.df)
summary(mod9)

mod9a<-phyloglm(pro~pol*heigh_height,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 30, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod9a)
