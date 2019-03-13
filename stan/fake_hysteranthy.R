
####March 4 2019/ Dan notes, the model is running but not really returning the proper parameters.

##wait, I now its returning the propper paramenters, wih a ton of divergetn transitions

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()
### fake hysteranthy data

#library(devtools)
#install_github("rmcelreath/rethinking")
library(rstan)
setwd("~/Documents/git/proterant/stan")


set.seed(613)
# Params
alpha<-1.5
b_pol<-3
b_flotime<--.5
b_minP<--.1


# x data
ndata <- 4000
pol<-rbinom(ndata,1,0.5)
flotime<-runif(ndata,4.5,8.5)
minP<-rnorm(ndata,8,2)

z <-alpha+b_pol*pol+b_flotime*flotime+b_minP*minP
# inverse logit -- this is the reverse of logit(p) in the model code
p <- 1/(1+exp(-z))
# now, can get bernoulli
y <- rbinom(ndata, 1, p)


dat <- as.data.frame(cbind(pol,flotime,minP,y))


#dat$pol.z<-NA
#dat$flo_time.z<-NA
#dat$


datalist<- with(dat, 
                  list(y=y,
                       pol=pol,
                       flotime=flotime,
                       minP=minP,
                  N = nrow(dat)
                             ))

berny<- stan('binary_stan_nophylo.stan', data = datalist,
                                        iter = 3000, warmup=2000) 
berny.sum <- summary(berny)$summary

berny.sum[c("alpha","b_pol","b_flotime","b_minP"),]
#This is from rethinking, but i should tray and do it in r stan
#flist <- alist(
 # y ~ dbinom(1, p),
#  logit(p) <-a + b_wind*pol + b_flo*flo_time+b_drought*min_p,
#  a <- dnorm(0, 3),
#  b_flo <- dnorm(0, 5),
#  b_wind <- dnorm(0, 5),
#  b_drought<-dnorm(0,5)
#)

#rbeta(200,1,1)
#mp.1 <- rethinking::map(flist, data = dat)

2.427140e-02
