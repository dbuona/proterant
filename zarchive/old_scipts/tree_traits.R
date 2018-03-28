###This makes a tree with all my variable. VERY TEDIOUS

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(tidyverse)
library(boot)
library(phylolm)

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")

####adjust red oak group
mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
#mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9  
mich.data$fruiting<-as.numeric(mich.data$fruiting)

####make variable
z <- as.factor(mich.data$pol); names(z) <- mich.tree$tip.label
zz<- as.factor(mich.data$pro); names(z) <- mich.tree$tip.label
zzz<-as.factor(mich.data$shade_bin); names(z) <- mich.tree$tip.label

mycol<-character(length(z))
mycol[mich.data$pro2==0]<-"green"
mycol[mich.data$pro2==1]<-"red"

mycol2<-character(length(zz))
mycol2[mich.data$pol==0]<-"yellow"
mycol2[mich.data$pol==1]<-"lightblue"

mycol3<-character(length(zzz))
mycol3[mich.data$shade_bin==0]<-"white"
mycol3[mich.data$shade_bin==1]<-"black"

fr <- mich.data$fruiting
names(fr) <- mich.tree$tip.label

fl <- mich.data$flo_time
names(fl) <- mich.tree$tip.label

H <- mich.data$heigh_height
names(H) <- mich.tree$tip.label

ourcol<-c("red","green","lightblue","yellow","white","black","blue")

par(mar=c(0,0,0,0))

plot(mich.tree,type = "fan", cex=0.4,,x.lim=c(0,300))
points(rep(210, length(mich.tree$tip.label)), 1:length(mich.tree$tip.label), pch=22, bg=mycol, cex=.7) 
text(200, length(mich.tree$tip.label)+1.5, "Hysteranthy",pos=4, cex=0.4)


points(rep(220, length(mich.tree$tip.label)), 1:length(mich.tree$tip.label), pch=22, bg=mycol2, cex=.7) 
text(213, length(mich.tree$tip.label)+1.5, "Syndrome",pos=4, cex=0.4)
points(rep(230, length(mich.tree$tip.label)), 1:length(mich.tree$tip.label), pch=22, bg=mycol3, cex=.7) 
text(225, length(mich.tree$tip.label)+1.5, "Shade",pos=4, cex=0.4)
segments(238, 1:nrow(mich.data), (238+fl), 1:nrow(mich.data), lwd=1, col="dark gray")
text(233, length(mich.tree$tip.label)+1.5, "Flower time", pos=4, cex=0.4)
segments(250, 1:nrow(mich.data), (250+fr), 1:nrow(mich.data), lwd=1, col="dark gray")
text(248, length(mich.tree$tip.label)+1.5, "Fruit time", pos=4, cex=0.4)
points(rep(270, nrow(mich.data)), 1:nrow(mich.data), pch=21, bg="blue", col="white", lwd=0.25, cex=H/12)
text(263, length(mich.tree$tip.label)+1.5, "Height", pos=4, cex=0.4)
legend(x = c(1, 60),y=c(110,130),legend=c("hysteranthous","non-hysteranthous","wind pollinated","insect pollinated","shade intolerant","shade tolerant","relative height"),col=ourcol,fill=ourcol, cex=0.5,
       box.lty=0)

####now do it for silvics

silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")

silv.data$fruiting<-NA
silv.data$fruiting<-silv.data$av_fruit_time
silv.data$fruiting[silv.data$fruiting==21]<-9

####make variable
z <- as.factor(silv.data$pol); names(z) <- silv.tree$tip.label
zz<- as.factor(silv.data$pro); names(z) <- silv.tree$tip.label
zzz<-as.factor(silv.data$shade_bin); names(z) <- silv.tree$tip.label

mycol<-character(length(z))
mycol[silv.data$pro==0]<-"green"
mycol[silv.data$pro==1]<-"red"

mycol2<-character(length(zz))
mycol2[silv.data$pol==0]<-"yellow"
mycol2[silv.data$pol==1]<-"lightblue"

mycol3<-character(length(zzz))
mycol3[silv.data$shade_bin==0]<-"white"
mycol3[silv.data$shade_bin==1]<-"black"

fr <- silv.data$fruiting
names(fr) <- silv.tree$tip.label

fl <- silv.data$flower_time
names(fl) <- silv.tree$tip.label

H <- silv.data$height
names(H) <- silv.tree$tip.label

ourcol<-c("red","green","lightblue","yellow","white","black","blue")

par(mar=c(1,0,0,1))

plot(silv.tree, cex=0.4,,x.lim=c(0,300))
points(rep(210, length(silv.tree$tip.label)), 1:length(silv.tree$tip.label), pch=22, bg=mycol, cex=.7) 
text(200, length(silv.tree$tip.label)+1.5, "Hysteranthy",pos=4, cex=0.4)
points(rep(220, length(silv.tree$tip.label)), 1:length(silv.tree$tip.label), pch=22, bg=mycol2, cex=.7) 
text(213, length(silv.tree$tip.label)+1.5, "Syndrome",pos=4, cex=0.4)
points(rep(230, length(silv.tree$tip.label)), 1:length(silv.tree$tip.label), pch=22, bg=mycol3, cex=.7) 
text(225, length(silv.tree$tip.label)+1.5, "Shade",pos=4, cex=0.4)
segments(238, 1:nrow(silv.data), (238+fl), 1:nrow(silv.data), lwd=1, col="dark gray")
text(233, length(silv.tree$tip.label)+1.5, "Flower time", pos=4, cex=0.4)
segments(250, 1:nrow(silv.data), (250+fr), 1:nrow(silv.data), lwd=1, col="dark gray")
text(248, length(silv.tree$tip.label)+1.5, "Fruit time", pos=4, cex=0.4)
points(rep(270, nrow(silv.data)), 1:nrow(silv.data), pch=21, bg="blue", col="white", lwd=0.25, cex=H/12)
text(263, length(silv.tree$tip.label)+1.5, "Height", pos=4, cex=0.4)
legend(x = c(1, 60),y=c(110,130),legend=c("hysteranthous","non-hysteranthous","wind pollinated","insect pollinated","shade intolerant","shade tolerant","relative height"),col=ourcol,fill=ourcol, cex=0.5,
       box.lty=0)
