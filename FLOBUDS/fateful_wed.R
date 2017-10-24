###fateful wednesday. Made by dan on Tuesday the 24 Oct. THis is what I should use when processing twigs.

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(splitstackshape)
library(reshape)
library(randomizr)

setwd("~/Documents/git/proterant/fLOBUDS")


#1organize all bags by species
dat<-read.csv("final_flobud_datapoints.csv")
dat$collected<-as.numeric(dat$collected)
sheet<-filter(dat,!is.na(collected))


sheet<-expandRows(sheet, "collected")

#try to obtain more cuttings for
###PRU VIR (10)
###ILE MUC 5
###ACE SAC 4
###Ace pen (15)  (or at least try to get to 42)

sheet$collect[sheet$name == "X"] <- #  ###if I cant get 48 of everything will need to sperate out the smaller treatrment
  
sheet<-expandRows(sheet, "collected") 

##THis randomly selects 48 of each species without replacement
by_cyl <- sheet %>% group_by(nomen)
exper<-sample_n(by_cyl, 48, replace = FALSE)


#assign to treatment
Z <- block_ra(block_var = pig$name, condition_names = c("WL0", "WS0", "WL1","WS1","CL0","CS0","CL1","CS1"))
exper$assignment<-Z
  
###check it
Ac<-filter(exper,nomen=="Acer_rubrum")
nrow(Ac) ###should be 48
table(Ac$assignment)

CP<-filter(exper,nomen=="Comptonia_peregrina")
nrow(Cc) ###should be 48
table(Cp$assignment)

###assign to beaker
X <- block_ra(block_var = exper$assignment,condition_names = 1:24)
exper$group<-X
table(exper$group)
exper$flask_id<-paste(exper$assignment,exper$group,sep="_")

write.csv() #final data sheet

tags<-dplyr::select(sheet name,flask_id)### this can be used to make name tags
write.csv() ### name tag output



