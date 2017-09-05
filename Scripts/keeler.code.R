###This is Dan's analysis of keeler dataset (1907) began 8/20.17

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
#anthy<-read.csv("keeler.csv", header = TRUE)

###Check names and clean
#library('Taxonstand')
#taxa <- paste(anthy$Genus,anthy$Species,sep=" ")
#taxa <- unique(taxa)
##matching to TPL 1.1
#clean_names <- TPL(taxa) # patience, patience
#clean_names$name<-paste(clean_names$New.Genus,clean_names$New.Species,sep="_")
#anthy$name<-clean_names$name
#clean_names<- filter(clean_names, Plant.Name.Index==TRUE )
#anthy<-filter(anthy, name %in% clean_names$name)
#write.csv(anthy, "cleaned_keeler.csv")

anthy<-read.csv("cleaned_keeler.csv")

anthy$name[anthy$name == "Polulus_grandidentata"] <- "Populus_grandidentata"
anthy$name[anthy$name == "Diospoyros_virginiana"] <- "Diospyros_virginiana"
anthy$name[anthy$name == "Ilex_monticola"] <- "Ilex_montana"
anthy$name[anthy$name == "Craetaegus_crus-gali"] <- "Crataegus_crus-gali"
anthy$name[anthy$name == "Craetaegus_coccinea"] <- "Crataegus_coccinea"
anthy$name[anthy$name == "Craetaegus_mollis"] <- "Crataegus_mollis"
anthy$name[anthy$name == "Craetaegus_punctata"] <- "Crataegus_punctata"
anthy$name[anthy$name == "Amelachier_canadensis"] <- "Amelanchier_canadensis"
anthy$name[anthy$name == "Quecus_coccinea"] <- "Quercus_coccinea"

##### Prune
anthy<-filter(anthy, !is.na(flo_time))

#write.csv(anthy,"cleaned_keeler.csv",row.names=FALSE)
names.intree<-treee$tip.label

# list of my species myspecies
namelist<-unique(anthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
#plot(pruned.by.anthy)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #60

addins<-setdiff(namelist,mytree.names)

###make ulta metric
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)

mytree.names<-pruned.by.anthy$tip.label


###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelistK<-final.df$name
namelistK==mytree.names
final.df$name== mytree.names

###data wrangle
final.df["pro"]<-NA
final.df$pro[final.df$Phen.sequence == "pro"] <- 1
final.df$pro[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$pro[final.df$Phen.sequence== "syn"] <- 0
final.df$pro[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$pro[final.df$Phen.sequence== "ser"] <- 0 
final.df$pro[final.df$Phen.sequence== "hyst"] <- 0

final.df["pro2"]<-NA
final.df$pro2[final.df$Phen.sequence == "pro"] <- 1
final.df$pro2[final.df$Phen.sequence == "pro/syn"] <- 1
final.df$pro2[final.df$Phen.sequence== "syn"] <- 1
final.df$pro2[final.df$Phen.sequence== "syn/ser"] <- 0
final.df$pro2[final.df$Phen.sequence== "ser"] <- 0 
final.df$pro2[final.df$Phen.sequence== "hyst"] <- 0


final.df["pol"]<-0 ### note, "ambo is coded as 0 right now"
final.df$pol[final.df$Pollination == "ambo"] <- 0
final.df$pol[final.df$Pollination == "insect"] <- 0
final.df$pol[final.df$Pollination == "wind"] <- 1

###convert height to meters
final.df$heigh_height.ft<-final.df$heigh_height.ft*.3048



##make it just shrub vs. tree
final.df$class2<-NA
final.df<- within(final.df, class2[heigh_height.ft<=15]<-0)
final.df<- within(final.df, class2[heigh_height.ft>15]<-1)

###reduce number variables for columns for flower class
final.df<- within(final.df, flower_class[flower_class=="polygamo-monecious"]<-"monoecious")
final.df<- within(final.df, flower_class[flower_class=="polygamo-dioecious"]<-"dioecious")
final.df<- within(final.df, flower_class[flower_class=="polygamous"]<-"perfect")

#dummy variable flower class
final.df$flo_type<-NA
final.df<- within(final.df, flo_type[flower_class=="perfect"]<-0)
final.df<- within(final.df, flo_type[flower_class=="monoecious"]<-1)
final.df<- within(final.df, flo_type[flower_class=="monecious"]<-1)
final.df<- within(final.df, flo_type[flower_class=="dioecious"]<-1)


### deal with av_fruit_time
hist(final.df$av_fruit_time)
fruits<-final.df$av_fruit_time

mean(final.df$av_fruit_time)
median(final.df$av_fruit_time)

final.df$fruit_bin<-NA
final.df<- within(final.df, fruit_bin[av_fruit_time<=8.5]<-0)
final.df<- within(final.df, fruit_bin[av_fruit_time>8.5]<-1)

final.df$pro<-as.integer(final.df$pro)
final.df$pro2<-as.integer(final.df$pro2)
final.df$pol<-as.integer(final.df$pol)
final.df$class2<-as.integer(final.df$class2)
final.df$shade_bin<-as.integer(final.df$shade_bin)
final.df$fruit_bin<-as.integer(final.df$fruit_bin)
final.df$flo_type<-as.integer(final.df$flo_type)


pruned.by.anthy$node.label<-NULL

write.csv(final.df,"keeler_cleaned.csv",row.names=FALSE)
write.tree(pruned.by.anthy,"pruned_keeler.tre")
goober<-na.omit(final.df)

##### Models


####full model 
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")
full.mod<-glm(pro~pol+flo_time,family = binomial(link="logit"),data=final.df)
summary(full.mod)



full.modA<-phyloglm(pro~pol+flo_time,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)
summary(full.modA)
