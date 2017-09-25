###This script creates a unique id for each twig that will be cult
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
table(dd$classification)

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
X <- block_ra(block_var = dd$assignment, num_arms = 24)
dd$group<-X


###test itWS0<-filter(dd,assignment=="WS0")
table(WS0$group)
######seems to work but wont be able to tell until there are only 48 individuals per treament.

