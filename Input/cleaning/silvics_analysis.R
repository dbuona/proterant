###Dan checks results of main analysis with data from USFS silvics. Main eefects only
##23 Aug 2017
#updated 27 sept 2017
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
anthy<-read.csv("..//Data/silvics_data.csv", header = TRUE)
#anthy<-filter(anthy, !is.na(av_fruit_time))
anthy<-na.omit(anthy)

###Prune tree
#list of species in tree
names.intree<-treee$tip.label

#dataformat it like Zanne
anthy$name<-paste(anthy$genus,anthy$species,sep="_")

#fix names
anthy$name[anthy$name == "Fraxinus_pensylvanica"] <- "Fraxinus_pennsylvanica"
anthy$name[anthy$name=="Ailthanthus_altissima"]<-"Ailanthus_altissima"
# list of my species myspecies
namelist<-unique(anthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)

mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #70 species include

addins<-setdiff(namelist,mytree.names) #14didn't make it

###make ultrametric
is.ultrametric(pruned.by.anthy)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)

species<-addins
for(i in 1:length(species)) pruned.by.anthy<-add.species.to.genus(pruned.by.anthy,species[i],
                                                                  where="root")
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #82 species include
setdiff(namelist,mytree.names) #2


##add them back
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Acer_barbatum",genus=NULL,where="random")
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Populus_heterophylla",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Ulmus_thomasii",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Salix_nigra",genus=NULL,where="random")
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Quercus_phellos",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Quercus_nuttallii",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Carya_laciniosa",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Carya_aquatica",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Acer_nigrum",genus=NULL,where="random") added
#pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Aesculus_octandra",genus=NULL,where="random") added

###make ultrametric (using mean path length smoothing, could also try penalized maximum likelihood with chronos())


mytree.names<-pruned.by.anthy$tip.label

###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelist3<-final.df$name
namelist3==mytree.names
final.df$name== mytree.names


setdiff(mytree.names,namelist3)

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
mean(final.df$av_fruit_time) #10
median(final.df$av_fruit_time) #9.5

final.df$fruit_bin<-NA
final.df<- within(final.df, fruit_bin[av_fruit_time<=8.5]<-0)
final.df<- within(final.df, fruit_bin[av_fruit_time>8.5]<-1)

final.df<- within(final.df, shade[shade=="very_tolerant"]<-"tolerant")
final.df<- within(final.df, shade[shade=="very_intolerant"]<-"intolerant")
final.df<- within(final.df, shade[shade=="medium"]<-"medium_tolerant")
unique(final.df$shade)

final.df$shade_bin<-NA
final.df<- within(final.df, shade_bin[shade=="medium_tolerant"]<-1)
final.df<- within(final.df, shade_bin[shade=="tolerant"]<-1)
final.df<- within(final.df, shade_bin[shade=="intolerant"]<-0)



##check again
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names


setdiff(mytree.names,namelist2)

pruned.by.anthy$node.label<-NULL

###output data sheets

write.csv(final.df, "silv_data_full.csv", row.names=FALSE)

write.tree(pruned.by.anthy,"pruned_silvics.tre")


########phylo signal###############
d<-comparative.data(pruned.by.anthy,final.df,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloD <- phylo.d(d, binvar=pro)##D=0.3156129

## models

final.df<-dplyr::select(final.df,name,pro,pol,av_fruit_time,flower_time,fruit_bin)
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")

silv.mod<-glm(pro~pol+flower_time+av_fruit_time,family = binomial(link="logit"),data=final.df)
summary(silv.mod)

silv.modcont<-phyloglm(pro~flower_time+pol+av_fruit_time,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)
summary(silv.modcont)

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

#####################################red oak group could be pulling things weird TRY TO change it back to annual
final.df$av_fruit_time[final.df$av_fruit_time == 22.0] <- 10.0
final.df$av_fruit_time[final.df$av_fruit_time == 21.0] <- 9.0


silv.mod<-glm(pro~pol+av_fruit_time,family = binomial(link="logit"),data=final.df)
summary(silv.mod)

silv.modcont<-phyloglm(pro~pol+flower_time+av_fruit_time,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                       start.beta=NULL, start.alpha=NULL,
                       boot=10,full.matrix = TRUE)
summary(silv.modcont)




#####
#final.df<-filter(final.df, av_fruit_time<= 14)
final.df<-rownames_to_column(final.df, "name")
inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

modelsilv <- brm(pro~pol+flower_time+av_fruit_time+(1|name), data = final.df, 
             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
             prior = c(prior(normal(0, 5), "b"),
                       prior(normal(0, 5), "Intercept"),
                       prior(student_t(3, 0, 5), "sd"))) 
summary(modelsilv)

############################### removing NAs
d<-read.csv("silv_data_full.csv",header=TRUE)
d<-na.omit(d)

t<-read.tree("pruned_silvics.tre")
namelist<-unique(d$name)
names.intree<-t$tip.label
##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(t,to.prune)

mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #55 species include

addins<-setdiff(namelist,mytree.names) #10didn't make it

###make ultrametric
is.ultrametric(pruned.by.anthy)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)
write.tree(t,"pruned_silv.tre")
