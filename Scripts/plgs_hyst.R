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

###now make pollination syndrom discrete
final.df["pol"]<-0
final.df$Pollination[final.df$Pollination == "wind"] <- "3wind"
final.df$Pollination[final.df$Pollination == "ambo"] <- "2ambo"
final.df$Pollination[final.df$Pollination == "insect"] <- "1insect"
final.df<- within(final.df, pol[is.na(Pollination)]<-1)
###now height
final.df$class<-NA
final.df<- within(final.df, class[heigh_height<10]<-"1shrub")
final.df<- within(final.df, class[heigh_height>=10 & heigh_height <20]<-"2small_tree")
final.df<- within(final.df, class[heigh_height>=20]<-"3large_tree")


###make factors for bianary traits
final.df$pro<-as.integer(final.df$pro)
final.df$pol<-as.integer(final.df$pol)
#final.df$class<-as.integer(final.df$class)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
head(final.df)
name.check(pruned.by.anthy, final.df)

#model
#pro <- final.df[, "pro"]
#pol<- final.df[, "pol"]
#height<-final.df[,"heigh_height"]

#plot variables
ggplot(final.df, aes(pro))+geom_bar()
ggplot(final.df,aes(pol))+geom_bar()
ggplot(final.df,aes(heigh_height))+geom_density()
ggplot(final.df,aes(class))+geom_bar()


library("phylolm")
##cleaning
final.df<- within(final.df, Pollination[is.na(Pollination)]<-"2ambo")
final.df<- within(final.df, Pollination[Pollination=="wnd"]<-"3wind")

##Dummy variable for pollination
final.df["pol"]<-NA
final.df$pol[final.df$Pollination == "3wind"] <- "3"
final.df$pol[final.df$Pollination == "2ambo"] <- "2"
final.df$pol[final.df$Pollination == "1insect"] <- "1"

###models: for 1 predictor: pollination syndrone
##How to I interprest this model with catagorical data? DO I need a dummy varaible? possibly helkpful https://onlinecourses.science.psu.edu/stat504/node/225
mod1<-glm(pro~Pollination-1,family = binomial(link="logit"),data=final.df)

summary(mod1)
library(boot)
inv.logit(-1.7492)
inv.logit(-0.1823)
inv.logit(-0.1911)
0.1481481 +0.4545508+ 0.4523699

mod2<-phyloglm(pro~as.factor-1,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
              start.beta=NULL, start.alpha=NULL,
              boot = 0, full.matrix = TRUE)

summary(mod2)
inv.logit(-1.65410)
inv.logit(0.30585)
inv.logit(-0.93966)
## check crawley ### height first
mod3<-phyloglm(pro~Pollination+class,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 20, log.alpha.bound = 4,
               start.beta=NULL, start.alpha=NULL,
               boot = 0, full.matrix = TRUE)
summary(mod3)

###interpret with divid by 4 rule?

