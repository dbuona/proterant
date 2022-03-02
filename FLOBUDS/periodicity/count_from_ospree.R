## Started 2 March 2019 ##
## By Lizzie ## and Dan

# Start your own R file to make a f(x) to count interactive experiments, ideally for all possible intxns across photoperiod, chilling, forcing, and field sample date (maybe provenance) treatments! #

## See also: https://github.com/lizzieinvancouver/ospree/issues/235

# housekeeping
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
library(dplyr)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses") 
} else setwd("~/Documents/git/ospree/analyses")

source("misc/getfielddates.R") # f(x) counts up field sample dates separated by a number of days you specify
source("misc/gettreatdists.R") # f(x) counts up treatment interactions, and more!


######Dan's summaries
#94 datasets
#152 experiments
#31 experiments manipulate photoperiod line 67
#8 of those have forcing periodicity suggesting 8/31 or ~25% of all photoperiod studies may have the issue

# or 40 experiments manipulate photoperiod based on line 122
# and 17 have the issue  42%

#107 experiments manipulate photoperiod and force (why is this number bigger than above) line 126
##I think I can answer that, line 126 mean 107 exp manipulated either photoperiod or forcing, even if the other was constatn

#15 experiments have interactions of photoperiod and forcing



###################
# All OSPREE data #
###################

dat <- read.csv("output/ospree_clean.csv", header = TRUE)
dat <- dat[dat$woody=="yes",]
dat$fieldsample.date <- as.Date(dat$fieldsample.date, format="%d-%b-%Y")
dat$doy <- format(dat$fieldsample.date, "%j")
########################################################################
####Q1### How many studies are in this version of ospree#################
#########################################################################
unique(dat$datasetID)A1a ### 94 datasets in this version of ospeee
nrow(dat %>% group_by(datasetID,study) %>% count()) ###A1b 152 experiments



# Get the number of field sampling dates that are 14 or more weeks apart, first for each datasetIDx study ...
ddatefx.all <- subset(dat, select=c("datasetID", "study", "fieldsample.date"))
ddatefx <- ddatefx.all[!duplicated(ddatefx.all), ]
ddatefx$datasetIDstudy <- paste(ddatefx$datasetID, ddatefx$study)

## Change main DF so only relevant dates are included ... not pretty, but should work
# First get the unique dates in a df
dates2weeks.count <- countfieldsample(ddatefx, 14)
uniquedates.df <- fieldsample.getuniquedates(ddatefx, 14)
uniquedates.df$selectcolumn <- paste(uniquedates.df$datasetIDstudy, uniquedates.df$date)
# Now subset to sane # of columnns
datsm <- subset(dat, select=c("datasetID", "study", "genus", "species", "forcetemp","forcetemp_night", "photoperiod_day", 
                              "fieldsample.date", "chilltemp", "chillphotoperiod", "chilldays"))
datsm$study
head(datsm)
 ### DAN sould source the cleaning code countinrtsvyion from count interavtionds source("limitingcues/source/countintxns_cleanosp.R")
source("~/Documents/git/proterant/FLOBUDS/periodicity/countintxns_cleanosp.R")
## Okay, formatting to look at intxns
datsm$force <- as.numeric(datsm$forcetemp)
datsm$forcenight<-as.numeric(datsm$forcetemp_night)
datsm$photo <- as.numeric(datsm$photoperiod_day)

datsm.noNA <- subset(datsm, is.na(force)==FALSE & is.na(photo)==FALSE & is.na(forcenight)==FALSE)


########################################################################
####Q2### How many studies manipulate photoperiod#################
#########################################################################

photodats<-datsm.noNA %>%                    
   group_by(datasetID,study) %>%         
  summarise(photoperiods = n_distinct(photo))


photodats<-filter(photodats,photoperiods>=2)
nrow(photodats) #A2 51 studies manipulate photoperiod (ie multiple treatments of photoperiod)


########################################################################
####Q2a### How many studies manipulate photoperiod and many covary it with thermo period #################
#########################################################################

photodots<-left_join(photodats,datsm.noNA)

photodots$thermop<-ifelse(photodots$forcetemp==photodots$forcetemp_night,"no","yes")

thermoprd<-dplyr::select(photodots,datasetID,study,thermop)
thermoprd<-distinct(thermoprd)
table(thermoprd$thermop) ##8 covary
22/51 ########A2a  have the issue




###


########################################################################
####Q3a### How many studies interactively manipulate forcing and photoperiod #################
#########################################################################
#lizzie's way
osp.fp <- get.treatdists(datsm.noNA, "photo", "force")
osp.fpintxn <- subset(osp.fp, intxn>=2)
osp.fpintxn[order(osp.fpintxn$datasetID),]

nrow(osp.fpintxn) #A 3a 18



####side bar#####
lookatunique <- get.uniquetreats(datsm.noNA, "photo", "force")

# Now (not pretty part) we'll take all NA dates ...
datsm$selectcolumn <- paste(datsm$datasetID, datsm$study, datsm$fieldsample.date)
datsm14d <- datsm[which(datsm$selectcolumn %in% uniquedates.df$selectcolumn),]

dim(datsm)
dim(datsm14d) # not such a huge loss of rows

# Check we don't lose any datasets!
setdiff(unique(paste(datsm$datasetID, datsm$study)), unique(paste(datsm14d$datasetID, datsm14d$study)))

datsm14d.noNA <- subset(datsm14d, is.na(force)==FALSE & is.na(photo)==FALSE)

photodats2<-datsm14d.noNA %>%                    
  group_by(datasetID,study) %>%         
  summarise(photoperiods = n_distinct(photo))

photodats2<-filter(photodats2,photoperiods>=2)
nrow(photodats2) #50

photodots2<-left_join(photodats2,datsm14d.noNA)

photodots2$thermop<-ifelse(photodots2$forcetemp==photodots2$forcetemp_night,"no","yes")

thermoprd2<-dplyr::select(photodots2,datasetID,study,thermop)
thermoprd2<-distinct(thermoprd2)
table(thermoprd2$thermop) ###### same answer is with other data






#osp14d.fp <- get.treatdists(datsm14d.noNA, "photo", "force")
#nrow(osp14d.fp) 
#unique(osp14d.fp$datasetID) 

#osp14d.fpintxn <- subset(osp14d.fp, intxn>=2) # 18experiments from 
#osp14d.fpintxn[order(osp14d.fpintxn$datasetID),]
#nrow(osp14d.fpintxn) ## 18 studies have interactions
#unique(osp14d.fpintxn$datasetID) #from 13 studies

#####ignore above becuse its the same answer as osp.fpin

#####END SIDE BAR################

###now indentfy which of the above might have periodicity issues
thermop<-dplyr::filter(datsm.noNA, datasetID %in%unique(osp.fpintxn$datasetID))

moreforcinginfo <- get.treatdists.daynight(thermop, "forcetemp", "forcetemp_night")


########################################################################
####Q3b### How many studies interactively manipulate forcing and photoperiod and have periodicity problems #################
#########################################################################
forcingvaried <- subset(moreforcinginfo, treatinfo!="forcing does not vary")
studiesinclconstantforce <- subset(forcingvaried, numconstantforce>0) # some studies have both
studiesinclforceperiodicity <- subset(forcingvaried, numdiffforce>0) 
nrow(studiesinclforceperiodicity) #10 out of 18 might have this issue
unique(studiesinclforceperiodicity$datasetID) #6

