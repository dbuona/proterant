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



###format the data in the same order as the tree
final.df<-anthy[match(mytree.names, anthy$name),]
namelist2<-final.df$name
namelist2==mytree.names
final.df$name== mytree.names



####add comlumns for analysis
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

final.df["pol"]<-0 ### note, "ambo is coded as 0 right now"
final.df$pol[final.df$Pollination == "insect"] <- 0
final.df$pol[final.df$Pollination == "wind"] <- 1


final.df$pro<-as.integer(final.df$pro)
final.df$pol<-as.integer(final.df$pol)
#final.df$class<-as.integer(final.df$class)

##make height catagorical
final.df$class<-NA
final.df<- within(final.df, class[heigh_height<10]<-0)
final.df<- within(final.df, class[heigh_height>=10 & heigh_height <20]<-1)
final.df<- within(final.df, class[heigh_height>=20]<-2)

##make it just shrub vs. tree
final.df$class2<-NA
final.df<- within(final.df, class2[heigh_height<=15]<-0)
final.df<- within(final.df, class2[heigh_height>15]<-1)


#rescale height to be in 10 meter increments
height10<-final.df$heigh_height/10

stop("this is not an error, just tells you we've finished the code")





