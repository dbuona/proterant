###This script creates a unique id for each twig that will be cult. ignore this and use fateful wednesday instead.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(splitstackshape)
library(reshape)
library(randomizr)

setwd("~/Documents/git/proterant/fLOBUDS")


d<-read.csv("flobud_datapoints.csv")

d$classification<-paste(d$Genus,d$species,sep="_")
table(d$classification)
#8 treatments, 6 cuttings per treatment within species= 48 total cuttings

###calculate how many cutting per indivudal (rounded up to nearest whole number)
cuttingper<-as.data.frame(table(d$classification))
cuttingper$cutting_number<-48/cuttingper$Freq
cuttingper<-filter(cuttingper, Freq>6)

d$cut<-NA
###number per idvidual/species +1 serving as a dummy varaible to be dropped
d$cut[d$classification=="Acer_pensylvanicum"]<-6
d$cut[d$classification=="Acer_rubrum"]<-5
d$cut[d$classification=="Acer_saccharum"]<-7
d$cut[d$classification=="Betula_allagheniensis"]<-6
d$cut[d$classification=="Comptonia_peregrina"]<-7
d$cut[d$classification=="Corylus_cornuta"]<-6
d$cut[d$classification=="Fagus_grandifolia"]<-6
d$cut[d$classification=="Ilex_mucronata"]<-8
d$cut[d$classification=="Ilex_verticilata"]<-5
d$cut[d$classification=="Prunus_pensylvanica"]<-6
d$cut[d$classification=="Prunus_virginiana"]<-8
d$cut[d$classification=="Vaccinium_corymbosum"]<-6
d$cut[d$classification=="Viburnum_acerifolium"]<-6


d<-filter(d, !is.na(cut)) ##3getr rid of other species

###duplicate the appropriate # of times so there willbe a row per idividual
d<-expandRows(d, "cut")

###give them a unique identifyer
dplus<-mutate(d,ind_id = make.unique((name)))

##drop the dummer variable so each unique id ends in a number
grepl("^.+(1|2|3|4|5|6|7|8)$",dplus$ind_id)
dd<-subset(dplus,grepl("^.+(1|2|3|4|5|6|7|8)$",ind_id))
dd<-dd[,c(13,1,2,3,4,14,5,6,7,8,9,10,11)]
dd$collected<-""

#### treatment assignment, can also use treatment_assign.R
projected<-as.data.frame(table(dd$classification))

combnum<-right_join(projected,cuttingper,by="Var1")
#write.csv(combnum,"cutting_perspecies_guide.csv",row.names = FALSE)
#write.csv(dd,"pre_sample_datasheet.csv", row.names = FALSE)

###This works,
Z <- block_ra(block_var = dd$species, condition_names = c("WL0", "WS0", "WL1","WS1","CL0","CS0","CL1","CS1"))
dd$assignment<-Z

###check it###############################
Ac<-filter(dd,species=="rubrum")
table(Ac$assignment)
cp<-filter(dd,species=="peregrina")
table(cp$assignment)
cc<-filter(dd,species=="cornuta")
table(cc$assignment)
#######################################seems good

###assign to beakers
X <- block_ra(block_var = dd$assignment,condition_names = 1:24)
dd$group<-X
table(dd$group)

###test it
WS0<-filter(dd,assignment=="WS0")
table(WS0$group)
######seems to work but wont be able to tell until there are only 48 individuals per treament.

table(dd$block_var)

dd$block_var <- with(dd, paste(assignment,group, sep = "_"))

dd<-group_by(dd,block_var)

label<-dplyr::select(dd,ind_id,block_var)


###assign to location in chamber
####code for placing them in chamber position##############3




position<-dd[!duplicated(dd$block_var),]

X <- block_ra(block_var = position$assignment,condition_names = 1:24)
position$placement<-X
position<-dplyr::select(position, assignment,group,block_var,placement)



WL0<-filter(position, assignment=="WL0")
table(WL0$block_var)
table(WL0$placement)
