
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

# Params
a <- 0.5
b_flo <- -.05
b_drought <- -0.001
b_wind<- .2


# x data
ndata <- 200


flo_time <- runif(ndata, 4.5, 9)
min_p <- rnorm(ndata, mean = 8, 1)
pol<-rbinom(ndata,size = 1,prob = 0.5)
# linear model
z <- a + b_wind*pol + b_flo*flo_time+b_drought*min_p
# inverse logit -- this is the reverse of logit(p) in the model code
p <- 1/(1+exp(-z))
# now, can get bernoulli
y <- rbinom(ndata, 1, p)

dat <- as.data.frame(cbind(pol, flo_time,min_p, y))
#dat$pol.z<-NA
#dat$flo_time.z<-NA
#dat$


datalist<- with(dat, 
                  list(y=y, 
                  pol = pol,
                    flo=flo_time,
                     minp=min_p, 
                  N = nrow(dat)
                             ))

berny<- stan('binary_stan_nophylo.stan', data = datalist,
                                        iter = 3000, warmup=2000) 

summary(berny)
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
