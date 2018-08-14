###PEP analy sis, use pep to find trends in hysteranthy
rm(list=ls()) 
options(stringsAsFactors = FALSE)
getOption("device")
options(device="quartz") ###if at arb if not RStudioGD
sessionInfo()

setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)
library(ggthemes)
library(brms)
big.dat<-read.csv("PEP_OFFSET_FULL.csv")

###does latitude or altitude predict offset

priorz<-get_prior(offset~alt+lat, data=big.dat)
cline.brms<-brm(offset~alt+lat+(alt+lat|taxa),data=big.dat, cores=4)
               

simple.cline<-brm(offset~alt+lat,data=big.dat)
