rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/proterant/input")

####cant get stan code to run
runstan = TRUE
mapsandsummaries = TRUE # these are just summaries for mapping (they're slow but not terribly), note that right now they run for the 10+ year dataset
use20yrs = FALSE # Otherwise it uses all data with 10 or more years
shinystancheck = FALSE # If you want to look at output in shiny stan
ncores = 2 # how many cores to use, only applies when runstan=TRUE

## Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(rgdal)
library(ggplot2)
library(lubridate)
library(rstan)
library(arm)
library(shinystan)

aln<-read.csv("alnus_delta_hyst.csv",header=TRUE)


alnALL <- aln[order(aln$s_id),]


#################################
## Look at the data and format ##
#################################

## Let's figure out sites and stages data ...
alnagg <- aggregate(alnALL[("year")], alnALL[c("s_id", "flower","leaf","offset", "lat", "lon", "alt")],
                    FUN=length)




# PEP_ID seems unique

# Subset the data based on the above for now ...

alnuse20 <- aln


####################################
## Get mean at each site and plot ##
####################################

## summarizing data
if(mapsandsummaries){
  meanaln <-
    ddply(alnALL, c("s_id", "lon", "lat"), summarise,
          mean = mean(offset),
          mean.yr = mean(year),
          sd = sd(offset),
          sem = sd(offset)/sqrt(length(offset)))
}  
 
pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)

}


alnUSW <- pepnumber(aln, "s_id")

#################################################
## Fit hinge models for each species (in Stan) ##
#################################################


# now add hinge
aln<-alnUSW
aln$YEAR.hin <- aln$year
aln$YEAR.hin[which(aln$YEAR.hin<1980)] <- 1980
aln$YEAR.hin <- aln$YEAR.hin-1980
# Note that I tried centering and scaling doy and it did not speed up much. The 20 yrs model still said: 1000 transitions using 10 leapfrog steps per transition would take 215.36 seconds.


# Note to self: for betula lmer will fit random intercepts but not random slopes

N <- nrow(aln)
y <- aln$offset
J <- length(unique(aln$peporder))
sites <- aln$peporder
year <- aln$YEAR.hin
  # nVars <-1
  # Imat <- diag(1, nVars)
  
  fit.hinge.aln <- stan("stan/hinge_randslopesint.stan",
                        data=c("N","J","y","sites","year"), iter=2500, warmup=1500,
                        chains=4, cores=ncores)
  # control = list(adapt_delta = 0.95, max_treedepth = 15))

}  
save(fit.hinge.aln, file="stan/fit.hinge.aln.Rda")
summary(fit.hinge.aln)
launch_shinystan(fit.hinge.aln)
class(fit.hinge.aln)
list_of_draws <- extract(fit.hinge.aln)
print(names(list_of_draws))

fit_summary <- summary(fit.hinge.aln)
print(names(fit_summary))
mu_tau_summary <- summary(fit.hinge.aln, pars = c("mu_b","sigma_b"), probs = c(0.1, 0.9))$summary
print(mu_tau_summary)
plot(fit.hinge.aln, pars="mu_b", show_density = TRUE, ci_level = 0.95, fill_color = "purple")
posterior_predict(fit.hinge.aln)
library(brms)
fit.aln.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aln) 
