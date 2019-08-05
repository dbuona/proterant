rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")
library(ggplot2)

HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
weather<-read.csv("..//FLOBUDS/data/hf000-01-daily-m.csv",header = TRUE)
head (weather)
unique(HF$species)
Vc<-filter(HF,species=="VACO")
Vc<-tidyr::gather(Vc, phase,DOY,4:7)
Vc1<-Vc %>% filter(phase %in% c("fopn.jd","l75.jd"))
Vc2<-Vc %>%filter(phase %in% c("fopn.jd","bb.jd"))
ggplot(Vc1,aes(year,DOY))+stat_summary(aes(color=phase))+theme_bw()
ggplot(Vc2,aes(year,DOY))+stat_summary(aes(color=phase))+theme_bw()
## we predict
#1990 91 should be low chilling or forcing
#'93 98 and 2000 high chilling or force
#'#94-96 could be low f and c with long photoperoi
#'
#'       
install.packages("chillR")             
library("chillR")
https://cran.r-project.org/web/packages/chillR/vignettes/hourly_temperatures.html