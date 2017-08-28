###This will be the source code for calling my phylogenetic tree adn michigan trees data to be used in further analysis

#list of species in tree
names.intree<-treee$tip.label

#dataformat it like Zanne
anthy$name<-paste(anthy$Genus,anthy$Species,sep="_")

# list of my species myspecies
namelist<-unique(anthy$name)

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
pruned.by.anthy<-drop.tip(treee,to.prune)
#plot(pruned.by.anthy)

###what are the tip labels in pruned phylogeny?
mytree.names<-pruned.by.anthy$tip.label

intersect(namelist,mytree.names) #107 species include

setdiff(namelist,mytree.names) #30 species did not make it
###
###make ultrametric (using mean path length smoothing, could also try penalized maximum likelihood with chronos())
is.ultrametric(pruned.by.anthy)
help(chronoMPL)
pruned.by.anthy<-chronoMPL(pruned.by.anthy)
is.ultrametric(pruned.by.anthy)
#plot(pruned.by.anthy)
#adding species to tree
pruned.by.anthy


pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Prunus_nigra",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Populus_heterophylla",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Quercus_coccinea",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Amelanchier_laevis",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Salix_pedicellaris",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Salix_nigra",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Ulmus_thomasii",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Quercus_ellipsoidalis",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Amelanchier_sanguinea",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Crataegus_succulenta",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Salix_serissima",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Carya_lacinosa",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Salix_candida",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Malus_coronaria",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Quercus_muehlenbergii",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Acer_nigrum",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Crataegus_crus-galli",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Salix_humilis",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Vaccinium_myrtilloides",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Sorbus_americana",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Quercus_prinoides",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Alnus_viridis",genus=NULL,where="random")
pruned.by.anthy<-add.species.to.genus(pruned.by.anthy, "Crataegus_macrosperma",genus=NULL,where="random")

mytree.names<-pruned.by.anthy$tip.label


###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names


setdiff(mytree.names,namelist2)


####add comlumns for analysis
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

###now make pollination syndrom bianary
##change acer
final.df$Pollination[final.df$name == "Acer_spicatum"] <- "insect"
final.df$Pollination[final.df$name == "Acer_pensylvanicum"] <- "insect"
final.df$Pollination[final.df$name == "Acer_negundo"] <- "wind"
final.df$Pollination[final.df$name == "Acer_platanoides"] <- "wind"
final.df$Pollination[final.df$name == "Acer_rubrum"] <- "wind"
final.df$Pollination[final.df$name == "Acer_saccharinum"] <- "wind"
final.df$Pollination[final.df$name == "Acer_saccharum"] <- "wind"

final.df["pol"]<-0 ### note, "ambo is coded as 0 right now"
final.df$pol[final.df$Pollination == "insect"] <- 0
final.df$pol[final.df$Pollination == "wind"] <- 1



##make height catagorical
final.df$class<-NA
final.df<- within(final.df, class[heigh_height<10]<-"0shrub")
final.df<- within(final.df, class[heigh_height>=10 & heigh_height <20]<-"1smalltree")
final.df<- within(final.df, class[heigh_height>=20]<-"2talltree")

##make it just shrub vs. tree
final.df$class2<-NA
final.df<- within(final.df, class2[heigh_height<=15]<-0)
final.df<- within(final.df, class2[heigh_height>15]<-1)

###reduce number variables for columns for flower class
final.df<- within(final.df, flower_class[flower_class=="polygamo_monecious"]<-"monoecious")
final.df<- within(final.df, flower_class[flower_class=="polygamo_dioecious"]<-"dioecious")
final.df<- within(final.df, flower_class[flower_class=="polygamous"]<-"perfect")

#dummy variable flower class
final.df$flo_type<-NA
final.df<- within(final.df, flo_type[flower_class=="perfect"]<-0)
final.df<- within(final.df, flo_type[flower_class=="monoecious"]<-1)
final.df<- within(final.df, flo_type[flower_class=="dioecious"]<-1)

##try collapsing tolerance
final.df<- within(final.df, shade_tol[shade_tol=="very_tolerant"]<-"tolerant")
final.df<- within(final.df, shade_tol[shade_tol=="very_intolerant"]<-"intolerant")
unique(final.df$shade_tol)

####option for bianry tolerance
final.df$shade_bin<-NA
final.df<- within(final.df, shade_bin[shade_tol=="moderately_tolerant"]<-1)
final.df<- within(final.df, shade_bin[shade_tol=="tolerant"]<-1)
final.df<- within(final.df, shade_bin[shade_tol=="intolerant"]<-0)

### deal with av_fruit_time
hist(final.df$av_fruit_time)
fruits<-final.df$av_fruit_time
final.df<-filter(final.df, !is.na(av_fruit_time))
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

#final.df<-dplyr::select(final.df, name,pro,pro2,pol,class2,flo_type,fruit_bin,shade_bin)

stop("this is not an error, just tells you we've finished the code")





