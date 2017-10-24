### Datasheet build--ignore this script and use fateful wednesdayt instead.
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(splitstackshape)
library(reshape)
library(randomizr)

setwd("~/Documents/git/proterant/fLOBUDS")


d<-read.csv("post_data_build.csv", header=TRUE)
d<-gather(d,arb,name,A:L)
d$ind_num<-" "
d$unique<-paste(d$name,d$ind_num,sep="_")

Z <- block_ra(block_var = d$name, condition_names = c("WL0", "WS0", "WL1","WS1","CL0","CS0","CL1","CS1"))
d$assignment<-Z


##CHeck it
Ac<-filter(d,name=="ACE_RUB")
table(Ac$assignment)
cp<-filter(d,name=="COM_PER")
table(cp$assignment)
cc<-filter(d,name=="COR_COR")
table(cc$assignment)

###assign to beakers
X <- block_ra(block_var = d$assignment,condition_names = 1:24)
d$group<-X
table(d$group)

###test it
WS0<-filter(d,assignment=="WS0")
table(WS0$group)

d$flask_id<-paste(d$assignment,d$group,sep="_")

write.csv(d,"label_making.csv",row.names = FALSE)

dd<-select(d,name,flask_id)


#### how to decide which to use
dat<-read.csv("final_flobud_datapoints.csv")
dat$collected<-as.numeric(dat$collected)
sheet<-filter(dat,!is.na(collected))

###up ACESAC # by splitting ACE_SACME 19 in 4 more


sheet<-expandRows(sheet, "collected")

###now the sheet has every twig on it
