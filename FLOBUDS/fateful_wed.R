###fateful wednesday. Made by Dan on Tuesday the 24 Oct. THis is what I should use when processing twigs.

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

#try to obtain more cuttings for
###PRU VIR (10)
###ILE MUC 5
###ACE SAC 4
###Ace pen (15)  (or at least try to get to 42)

sheet$collected[sheet$name == "PRUVIR 4 HF DB"] <-9
sheet$collected[sheet$name == "ACEPEN19 HF DB"] <-6
sheet$collected[sheet$name == "ACEPEN18 HF DB"] <-8
sheet$collected[sheet$name == "ACEPEN14 HF DB"] <-14
sheet$collected[sheet$name == "PRUVIR 3 HF DB"] <-5
sheet$collected[sheet$name == "PRUVIR 5 HF DB"] <-9  
sheet$collected[sheet$name == "PRUVIR 7 HF DB"] <-9 
sheet$collected[sheet$name == "ACEPEN17 HF DB"] <-10
sheet$collected[sheet$name == "ACESAC19 HF DB"] <-10
sheet$collected[sheet$name == "ACESAC18 HF DB"] <-10
sheet$collected[sheet$name == "ACEPEN20 HF DB"] <-14
sheet$collected[sheet$name == "ACEPEN15 HF DB"] <-6

sheet$collected[sheet$name == "ILEMUC 12 HF DB"] <-10
sheet$collected[sheet$name == "ILEMUC 14 HF DB"] <-8
sheet$collected[sheet$name == "ILEMUC 13 HF DB"] <-9
sheet$collected[sheet$name == "ILEMUC 11 HF DB"]<-8


sheet<-expandRows(sheet, "collected") 

##THis randomly selects 48 of each species without replacement
by_cyl <- sheet %>% group_by(nomen)
set.seed(1800)
exper<-sample_n(by_cyl, 48, replace = FALSE)
useful<-table(exper$name)
useful<-as.data.frame(useful)
write.csv(useful,"twigstouseinexp.csv",row.names = FALSE)

####At this point, I will print the above made data sheet lab help goes and pulls the proper # of individuals from each bag, 
###defoliate and put in new bags

#assign to treatment
Z <- block_ra(block_var = exper$nomen, condition_names = c("WL0", "WS0", "WL1","WS1","CL0","CS0","CL1","CS1"))
exper$assignment<-Z
  
###check it
Ac<-filter(exper,nomen=="Acer_rubrum")
nrow(Ac) ###should be 48
table(Ac$assignment)

CP<-filter(exper,nomen=="Comptonia_peregrina")
nrow(CP) ###should be 48
table(CP$assignment)

###assign to beaker
X <- block_ra(block_var = exper$assignment,condition_names = 1:24)
exper$group<-X
table(exper$group)
exper$flask_id<-paste(exper$assignment,exper$group,sep="_")

write.csv(exper,"final.data.sheet.csv",row.names = FALSE) #final data sheet

tags<-dplyr::select(exper, name,flask_id)### this can be used to make name tags
write.csv(tags,"flask_labels.csv")### name tag output

### update to avoid con specifics
lapel<-read.csv("flask_labels.csv",header=TRUE)
lapel$good_flaskid<-NA
changes<-filter(lapel, !is.na(new_num))
changes$good_flaskid<-paste(changes$treat,changes$new_num,sep="_")
nochange<-filter(lapel, is.na(new_num))
nochange$good_flaskid<-nochange$flask_id
test<-rbind(changes,nochange)

uselabel<-dplyr::select(test,name,good_flaskid)

write.csv(uselabel,"final_labels.csv",row.names = FALSE)

###Fill flasks with DI water
### color 



