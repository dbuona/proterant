### reshuffling ### this was probably the least effecient way to do this task, 
#which I am not even sure why I needed to do it, but Lizzie and I will try to remember
# Any way Norvermber 13 2017


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
library(ggplot2)

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")

silv.tree<-read.tree("pruned_silvics.tre")
silv.data<-read.csv("silv_data_full.csv")


mich.data$av_fruit_time[mich.data$av_fruit_time=="persistant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="persitant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="unreported"]<-9    
mich.data$av_fruit_time<-as.numeric(mich.data$av_fruit_time)

###clean fruiting
#mich.data$fruiting[mich.data$fruiting==19]<-7

#####add a new column for a adjusting for red acorn time
mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9                                      
mich.data$fruiting<-as.numeric(mich.data$fruiting)

####Silvics cleaning
###fruiting
silv.data$fruiting<-NA
silv.data$fruiting<-silv.data$av_fruit_time
silv.data$fruiting[silv.data$fruiting==21]<-9

###functional hysteranthy
silv.data["pro2"]<-NA
silv.data$pro2[silv.data$silvic_phen_seq== "pro"] <- 1
silv.data$pro2[silv.data$silvic_phen_seq== "pro/syn"] <- 1
silv.data$pro2[silv.data$silvic_phen_seq== "syn"] <- 1
silv.data$pro2[silv.data$silvic_phen_seq== "syn/ser"] <- 0
silv.data$pro2[silv.data$silvic_phen_seq== "ser"] <- 0 
silv.data$pro2[silv.data$silvic_phen_seq== "hyst"] <- 0
silv.data$pro2[silv.data$name == "Quercus_laurifolia"] <- 1

#####Centering
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))


silv.data$height_cent<-(silv.data$height-mean(silv.data$height))/(2*sd(silv.data$height))
silv.data$fruit_cent<-(silv.data$fruiting-mean(silv.data$fruiting))/(2*sd(silv.data$fruiting))
silv.data$flo_cent<-(silv.data$flower_time-mean(silv.data$flower_time))/(2*sd(silv.data$flower_time))




mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<-  silv.data %>% remove_rownames %>% column_to_rownames(var="name")
###
min(mich.data$flo_time)
max(mich.data$flo_time)

####or run  hyst final analysis to here
nsp<-nrow(mich.data)

set.seed(202)
##original
mich5cent<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(mich5cent)

####shuffle 1
for (i in 1:nsp){
  mich.data$add1<-sample(-1:1,size = nsp,replace=TRUE)
  mich.data$shuf1<-mich.data$flo_time+mich.data$add1
  }
mich.data$shuf1_cent<-(mich.data$shuf1-mean(mich.data$shuf1))/(2*sd(mich.data$shuf1))

mich5shuf1<-phyloglm(pro~pol+height_cent+shuf1_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(mich5shuf1)

####shuffle 2
for (i in 1:nsp){
  mich.data$add2<-sample(-2:2,size = nsp,replace=TRUE)
  mich.data$shuf2<-mich.data$flo_time+mich.data$add2
}
mich.data$shuf2_cent<-(mich.data$shuf2-mean(mich.data$shuf2))/(2*sd(mich.data$shuf2))

mich5shuf2<-phyloglm(pro~pol+height_cent+shuf2_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=50,full.matrix = TRUE)
summary(mich5shuf2)
########shuffle 3
for (i in 1:nsp){
  mich.data$add3<-sample(-3:3,size = nsp,replace=TRUE)
  mich.data$shuf3<-mich.data$flo_time+mich.data$add3
}
mich.data$shuf3_cent<-(mich.data$shuf3-mean(mich.data$shuf3))/(2*sd(mich.data$shuf3))

mich5shuf3<-phyloglm(pro~pol+height_cent+shuf3_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=50,full.matrix = TRUE)
summary(mich5shuf3)
#####Now it gets more complicated because we dont want things going below zero or above 12
for (i in 1:nsp){
  ifelse(mich.data$flo_time<4,mich.data$add4<-sample(-3:3,size = nsp,replace=TRUE),mich.data$add4<-sample(-4:4,size = nsp,replace=TRUE))
  mich.data$shuf4<-mich.data$flo_time+mich.data$add4
}
mich.data$shuf4_cent<-(mich.data$shuf4-mean(mich.data$shuf4))/(2*sd(mich.data$shuf4))

mich5shuf4<-phyloglm(pro~pol+height_cent+shuf4_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=50,full.matrix = TRUE)
summary(mich5shuf4)
####Shuffle 5

for (i in 1:nsp){
  ifelse(mich.data$flo_time<5,mich.data$add5<-sample(-4:4,size = nsp,replace=TRUE),mich.data$add5<-sample(-5:5,size = nsp,replace=TRUE))
  mich.data$shuf5<-mich.data$flo_time+mich.data$add5
}
mich.data$shuf5_cent<-(mich.data$shuf5-mean(mich.data$shuf5))/(2*sd(mich.data$shuf5))

mich5shuf5<-phyloglm(pro~pol+height_cent+shuf5_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=50,full.matrix = TRUE)
summary(mich5shuf5)

####shuffle 6:
for (i in 1:nsp){
  ifelse(mich.data$flo_time<6,mich.data$add6<-sample(-5:5,size = nsp,replace=TRUE),mich.data$add6<-sample(-6:6,size = nsp,replace=TRUE))
  mich.data$shuf6<-mich.data$flo_time+mich.data$add6
}
mich.data$shuf6_cent<-(mich.data$shuf6-mean(mich.data$shuf6))/(2*sd(mich.data$shuf6))

mich5shuf6<-phyloglm(pro~pol+height_cent+shuf6_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=50,full.matrix = TRUE)
summary(mich5shuf6)


test<-dplyr::select(mich.data,flo_time,shuf1,shuf2,shuf3,shuf4,shuf5,shuf6)

#####plot thus
a<-as.data.frame(coef(mich5cent))
a<-rownames_to_column(a, "name")

b<-as.data.frame(coef(mich5shuf1))
b<-rownames_to_column(b, "name1")

c<-as.data.frame(coef(mich5shuf2))
c<-rownames_to_column(c, "name2")

d<-as.data.frame(coef(mich5shuf3))
d<-rownames_to_column(d, "name3")

e<-as.data.frame(coef(mich5shuf4))
e<-rownames_to_column(e, "name4")

f<-as.data.frame(coef(mich5shuf5))
f<-rownames_to_column(f, "name5")

g<-as.data.frame(coef(mich5shuf6))
g<-rownames_to_column(g, "name6")
H<-cbind(a,b)
I<-cbind(c,d)
J<-cbind(e,f)
K<-cbind(H,I)
L<-cbind(J,g)
plotme<-cbind(K,L)
plotme<-dplyr::select(plotme,-name1,-name2,-name3,-name4,-name5,-name6)
PLOT3<-gather(plotme,"model","coeficient",2:8)
ggplot(PLOT3,aes(model,coeficient,color=name,group=name))+geom_point()+geom_line(aes(color=name))

#####or
median(mich.data$flo_time)
mean(mich.data$flo_time)

earlymich<-filter(mich.data,flo_time<=5.5)
earlymich<-na.omit(earlymich)
###prune tree
mich.data<-rownames_to_column(mich.data, "name")
namelist<-earlymich$name
names.intree<-mich.tree$tip.label

##Prune the tree
to.prune<-which(!names.intree%in%namelist)
earlytree<-drop.tip(mich.tree,to.prune)

earlytree$tip.label==earlymich$name

#recenter
earlymich$height_cent<-(earlymich$heigh_height-mean(earlymich$heigh_height))/(2*sd(earlymich$heigh_height))
earlymich$fruit_cent<-(earlymich$fruiting-mean(earlymich$fruiting))/(2*sd(earlymich$fruiting))
earlymich$flo_cent<-(earlymich$flo_time-mean(earlymich$flo_time))/(2*sd(earlymich$flo_time))
earlymich$pol_cent<-(earlymich$pol-mean(earlymich$pol))/(2*sd(earlymich$pol))
earlymich$av_fruit_time_cent<-(earlymich$av_fruit_time-mean(earlymich$av_fruit_time))/(2*sd(earlymich$av_fruit_time))

earlymich<-  earlymich %>% remove_rownames %>% column_to_rownames(var="name")

earlycent<-phyloglm(pro~pol+height_cent+flo_cent+fruit_cent+shade_bin,earlymich, earlytree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                    start.beta=NULL, start.alpha=NULL,
                    boot=50,full.matrix = TRUE)
summary(earlycent)

###flo time still super important even when deal with just early

