## Dan want to find papers with different temperature
# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

# libraries
library(shinystan)
library(reshape2)
library(rstan)
library(rstanarm)
# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges") 
} else if (length(grep("ailene", getwd()))>0) {setwd("~/Documents/GitHub/ospree/analyses/ranges")
}else if(length(grep("Ignacio", getwd()))>0) { 
  setwd("~/GitHub/ospree/analyses/ranges") 
} else if(length(grep("catchamberlain", getwd()))>0) { 
  setwd("~/Documents/git/ospree/analyses/ranges") 
} else if(length(grep("danielbuonaiuto", getwd()))>0) { 
  setwd("~/Documents/git/ospree/analyses/ranges") 
}else setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")


######################################
# Flags to choose for bbstanleadin.R #
######################################

# Our flags for ranges, for now ... (see issue #379)
use.chillports = FALSE
use.zscore = TRUE
use.allspp = FALSE
use.multcuespp = FALSE
use.cropspp = FALSE
use.expramptypes.fp = FALSE
use.exptypes.fp = FALSE
use.expchillonly = FALSE


setwd("..//bb_analysis")
source("source/bbstanleadin.R")

colnames(bb.all)
bb.expphotoforce$forcetemp<-as.numeric(bb.expphotoforce$forcetemp)
bb.expphotoforce$forcetemp_night<-as.numeric(bb.expphotoforce$forcetemp_night)

photo<-bb.expphotoforce %>% group_by(datasetID) %>% distinct(photoperiod_day)%>% count()
force<-bb.expphotoforce %>% group_by(datasetID) %>% distinct(forcetemp)%>% count()
manip<-filter(photo,n>1)
manif<-filter(force,n==1)
intersect(manip$datasetID,manif$datasetID)

intersect(manip$datasetID,periodio$datasetID)


periodio<-bb.expphotoforce%>% filter(forcetemp!=forcetemp_night)

both<-periodio %>% group_by(datasetID) %>% distinct(photoperiod_day)%>% count()

unique(periodio$datasetID)
unique(bb.expphotoforce$datasetID)
