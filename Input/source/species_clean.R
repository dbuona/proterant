##began 21 Aug 2017 by Dan
###trying to check species names to make sure  Zanne tree isn't losing thing.

###after running I don't think anything helps.

install.packages("Taxonstand")
library('Taxonstand')

taxa <- paste(anthy$Genus,anthy$Species,sep=" ")
taxa <- unique(taxa)


##matching to TPL 1.1
clean_names <- TPL(taxa) # patience, patience

clean_names<- filter(clean_names, Plant.Name.Index==TRUE )
clean_names$name<-paste(clean_names$New.Genus,clean_names$New.Species,sep="_")


setdiff(anthy$name,clean_names$name)

goober<-filter(anthy, name %in% clean_names$name)

anthy<-goober


