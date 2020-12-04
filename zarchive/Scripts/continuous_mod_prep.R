

#########READ IN ALL DATA AND ASSOCIATED TREES##################



####code offset as binary



### This file is for prepping USDA traits for merging with hf
traits<-read.csv("HarvardForest/HF.species.trait.data.csv",header=TRUE)
traits<-dplyr::filter(traits,species!="QUAL")
###make things numeric
traits$pol<-ifelse(traits$pollination=="insect",0,1)
traits$flo_view<-ifelse(traits$fl_conspic=="Y",0,1)
traits$frost_free<-as.numeric(traits$frost_free)

#rescale predictors
#traits$height_cent<-(traits$height_mat-mean(traits$height_mat))/(2*sd(traits$height_mat))
#traits$cent_frost<-(traits$frost_free-mean(traits$frost_free))/(2*sd(traits$frost_free))
#traits$cent_minT<-(traits$min_temp-mean(traits$min_temp))/(2*sd(traits$min_temp))
#traits$cent_minP<-(traits$min_precip-mean(traits$min_precip))/(2*sd(traits$min_precip))
#traits$cent_roots<-(traits$root_depth_in-mean(traits$root_depth_in))/(2*sd(traits$root_depth_in))
#traits$cent_pol<-(traits$pol-mean(traits$pol))/(2*sd(traits$pol))
#traits$cent_flo_view<-(traits$flo_view-mean(traits$flo_view))/(2*sd(traits$flo_view))
#traits$seed_cent<-(traits$seed_pound-mean(traits$seed_pound,na.rm=TRUE))/(2*sd(traits$seed_pound,na.rm=TRUE))

###prep the tree
treee<-read.tree("HarvardForest/Vascular_Plants_rooted.dated.tre")
names.intree<-treee$tip.label
namelist<-unique(traits$name)

to.prune<-which(!names.intree%in%namelist)
HF.tree.pruned<-drop.tip(treee,to.prune)
mytree.names<-HF.tree.pruned$tip.label

#setdiff(namelist,mytree.names) ###Viburum
#intersect(namelist,mytree.names)

###make ultrametric
HF.tree.pruned<-chronoMPL(HF.tree.pruned)
is.ultrametric(HF.tree.pruned)

###add one species
HF.tree.pruned<-add.species.to.genus(HF.tree.pruned, "Viburnum_lantanoides",genus=NULL,where="root")

##matchthem ### not sure if I need to do this for brms
df<-traits[match(mytree.names, traits$name),]
namelist<-df$name
namelist==mytree.names
df$name== mytree.names
#### cool
HF.tree.pruned$node.label<-NULL


