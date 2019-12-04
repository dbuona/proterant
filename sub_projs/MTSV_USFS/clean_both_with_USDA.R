rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
mich.tree<-read.tree("datasheets_derived/MTSV_USFS/pruned_for_mich.tre")
mich.data<-read.csv("datasheets_derived/MTSV_USFS/mich_data_full_clean.csv")
USDA.dat.MTSV<-read.csv("..//Data/USDA_traitfor_MTSV.csv",header=TRUE)

##quick clean from before
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)

###make the tree work 
mich.tree$node.label<-NULL

###Prune tree to species with USDA data available
mich.data<-dplyr::left_join(mich.data,USDA.dat.MTSV,by="name")

######prune the tree for drought modeling
mich.data<-dplyr::filter(mich.data,!is.na(min._precip))
####prune tree to match reduced dataset
names.intree<-mich.tree$tip.label
namelist<-unique(mich.data$name)
to.prune<-which(!names.intree%in%namelist)
mich.tree.droughtprune<-drop.tip(mich.tree,to.prune)
mytree.names<-mich.tree.droughtprune$tip.label
### the smaller new tree is called mich.tree.droughtprune

#####clean USFS now
silv.tree<-read.tree("../sub_projs/MTSV_USFS/pruned_silvics.tre")
silv.data<-read.csv("../sub_projs/MTSV_USFS/silv_data_full.csv")
silv.USDA<-read.csv("silv.USDA.csv")
silv.tree$node.label<-NULL
silv.data<-left_join(silv.data,silv.USDA) ###maybe you should make a drought species list specifically for silvics to lose less species

####Silvics cleaning########
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

silv.data<-filter(silv.data,!is.na(min._precip)) ### get rid of nas

names.intree<-silv.tree$tip.label
namelist<-unique(silv.data$name)
to.prune<-which(!names.intree%in%namelist)
silv.tree.droughtprune<-drop.tip(silv.tree,to.prune)
mytree.names<-silv.tree.droughtprune$tip.label
setdiff(namelist,mytree.names) 
intersect(namelist,mytree.names)

write.csv(silv.data,"../sub_projs/MTSV_USFS/silvdata_final.csv")
write.tree(silv.tree.droughtprune,"../sub_projs/MTSV_USFS/silvtre_final.tre")

