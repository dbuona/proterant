###Dan checks results of main analysis with data from USFS silvics. Main eefects only
##23 Aug 2017
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
#https://academic-oup-com.ezp-prod1.hul.harvard.edu/sysbio/article-lookup/doi/10.1093/sysbio/syp074 Garland and Ives 2010

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("silvics_data.csv", header = TRUE)
anthy<-filter(anthy, !is.na(av_fruit_time))

###Prune tree
#list of species in tree
names.intree<-treee$tip.label

#dataformat it like Zanne
anthy$name<-paste(anthy$genus,anthy$species,sep="_")

# list of my species myspecies
namelist<-unique(anthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)

mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #53 species include

setdiff(namelist,mytree.names) #11 didn't make it

###make ultrametric (using mean path length smoothing, could also try penalized maximum likelihood with chronos())
is.ultrametric(pruned.by.anthy)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)
#plot(pruned.by.anthy)

mytree.names<-pruned.by.anthy$tip.label

###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names


setdiff(mytree.names,namelist2)

### add analysis columns

final.df["pro"]<-NA
final.df$pro[final.df$silvic_phen_seq== "pro"] <- 1
final.df$pro[final.df$silvic_phen_seq== "pro/syn"] <- 1
final.df$pro[final.df$silvic_phen_seq== "syn"] <- 0
final.df$pro[final.df$silvic_phen_seq== "syn/ser"] <- 0
final.df$pro[final.df$silvic_phen_seq== "ser"] <- 0 
final.df$pro[final.df$silvic_phen_seq== "hyst"] <- 0
final.df$pro[final.df$name == "Quercus_laurifolia"] <- 1

### update pollination
final.df$silvics_pollin[final.df$name == "Sassafras_albidum"] <- "insect"
final.df$silvics_pollin[final.df$name == "Juglans_nigra"] <- "wind"
final.df$silvics_pollin[final.df$name == "Celtis_laevigata"] <- "wind"
final.df$silvics_pollin[final.df$name == "Gleditsia_triacanthos"] <- "insect"
final.df$silvics_pollin[final.df$name == "Diospyros_virginiana"] <- "insect"
final.df$silvics_pollin[final.df$name == "Halesia_carolina"] <- "insect"
final.df$silvics_pollin[final.df$name == "Fraxinus_profunda"] <- "wind"
final.df$silvics_pollin[final.df$name == "Fraxinus_latifolia"] <- "wind"
final.df$silvics_pollin[final.df$name == "Fraxinus_nigra"] <- "wind"

final.df["pol"]<-0 
final.df$pol[final.df$silvics_pollin == "insect"] <- 0
final.df$pol[final.df$silvics_pollin == "wind"] <- 1

##fruit time
final.df$av_fruit_time[final.df$av_fruit_time == "late_spring/early_summer"] <- 5.5
final.df$av_fruit_time[final.df$av_fruit_time == "as late as midwinter"] <- 13
final.df$av_fruit_time[final.df$av_fruit_time == "fall/winter"] <- 11
final.df$av_fruit_time[final.df$av_fruit_time == "winter"] <- 12

final.df$av_fruit_time<-as.numeric(final.df$av_fruit_time)

##check again
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names


setdiff(mytree.names,namelist2)

pruned.by.anthy$node.label<-""

## models

#final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
final.df<-dplyr::select(final.df,name,pro,pol,av_fruit_time)

silv.mod<-glm(pro~pol+av_fruit_time,family = binomial(link="logit"),data=final.df)
summary(silv.mod)

full.modA<-phyloglm(pro~pol+av_fruit_time,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)

#Error in phyloglm(pro ~ pol + av_fruit_time, final.df, pruned.by.anthy,  : 
#the number of rows in the data does not match the number of tips in the tree.
#but yes it does
pruned.by.anthy$tip.label==final.df$name
pruned.by.anthy$tip.label!=final.df$name
#see?

summary(full.modA)

###let just see it in bayesian
library("brms")
library("MCMCglmm")

final.df<-rownames_to_column(final.df, "name")
inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

model <- brm(pro~pol+av_fruit_time+(1|name), data = final.df, 
             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=6000,
             prior = c(prior(normal(0, 5), "b"),
                       prior(normal(0, 5), "Intercept"),
                       prior(student_t(3, 0, 5), "sd"))) 
summary(model)
###fruit time is not significant (this is because possibly the range is too big)

### red oak group seem to be nullifying fruit time effect--try it without them

final.df<-filter(final.df, av_fruit_time<= 14)
inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

modelminusoak <- brm(pro~pol+av_fruit_time+(1|name), data = final.df, 
             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
             prior = c(prior(normal(0, 5), "b"),
                       prior(normal(0, 5), "Intercept"),
                       prior(student_t(3, 0, 5), "sd"))) 
summary(modelminusoak)
