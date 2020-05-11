rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_projs/")

library(ape)
library(phytools)
library(brms)
library(tibble)
library(ggstance)
library(ggplot2)
library("dplyr")
library("jpeg")
library("phylolm")
library(ggstance)


##read in the data
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
HFsubber<-read.csv("HarvardForest/HFdata4modeling.csv",header=TRUE)
HF.tree<-read.tree("HarvardForest/HFtree4modeling.tre")
###make fls measure
HF$phys.fls<-HF$bb.jd-HF$fbb.jd
HF$funct.fls<-HF$l75.jd-HF$fopn.jd
HF$inter.fls<-HF$bb.jd-HF$fopn.jd

###make catagorical FLS
HF$hyst.funct<-ifelse(HF$funct.fls>0,1,0)
HF$hyst.phys<-ifelse(HF$phys.fls>0,1,0)
HF$hyst.inter<-ifelse(HF$inter.fls>0,1,0)

### prune the tree
HF<-dplyr::filter(HF,species!=("QUAL")) ## quercus alba has no flowers
spforcontmods<-HFsubber$species ##subset of species good for this analysis

HF.data<-dplyr::filter(HF,species %in% c(spforcontmods))
HF.data<-dplyr::left_join(HF.data,HFsubber, by="species") ###This is the data for the continuous models
##zscore predictors for these models
HF.data$pol_cent<-(HF.data$pol-mean(HF.data$pol,na.rm=TRUE))/(2*sd(HF.data$pol,na.rm=TRUE))
HF.data$precip_cent<-(HF.data$min_precip-mean(HF.data$min_precip))/(2*sd(HF.data$min_precip))
HF.data$flo_cent<-(HF.data$fopn.jd-mean(HF.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.data$fopn.jd,na.rm=TRUE))

HF.data$flo_cent.neg<--(HF.data$flo_cent)
HF.data$precip_cent.neg<--(HF.data$precip_cent)

unique(HF.data$name)


freez1<-read.csv("HarvardForest/Choatetal2012.csv")
freez1$MPDQ..mean.precipitation.of.the.driest.quarter.<-as.numeric(freez1$MPDQ..mean.precipitation.of.the.driest.quarter.)
freez1$AI..aridity.index.<-as.numeric(freez1$AI..aridity.index.)


freez2<-read.csv("HarvardForest/Bartlettetal2012.csv")
freez1$name <- sub(" ", "_", freez1$Species)
freez2$name <- sub(" ", "_", freez2$Species)


try<-read.csv("HarvardForest/try.csv")
try$name <- sub(" ", "_", try$name)


MTSV<-read.csv("MTSV_USFS/michdata_final.csv")
USFS<-read.csv("MTSV_USFS/silvdata_final.csv")

droughtmetrics<-data.frame("Overlapping species entries"=c("Original","Choate 2012","Bartlett 2012","TRY"),
MTSV=c(length(unique(MTSV$name)),
length(intersect(MTSV$name,freez1$name)),
length(intersect(MTSV$name,freez2$name)),
length(intersect(MTSV$name,try$name))),
USFS=c(length(unique(USFS$name)),
       length(intersect(USFS$name,freez1$name)),
       length(intersect(USFS$name,freez2$name)),
       length(intersect(USFS$name,try$name))),
HF=c(length(unique(HF.data$name)),
          length(intersect(HF.data$name,freez1$name)),
          length(intersect(HF.data$name,freez2$name)),
          length(intersect(HF.data$name,try$name)))
)



###but do they correlate?
colnames(freez1)

freezish<-dplyr::select(freez1, name,ψ50,AI..aridity.index.,ψ50.safety.margin)
freezish$ψ50<-as.numeric(freezish$ψ50)
freezish$AI..aridity.index.<-as.numeric(freezish$AI..aridity.index.)
freezish$ψψ50.safety.margin<-as.numeric(freezish$ψ50.safety.margin)





#### so not alot of sepceis coverage.

##############now correlations between flowering time and other predictors
addmich<-read.csv("mich_tree_additions.csv")
addmich$flo_loc_num<-ifelse(addmich$flower_local=="lateral",1,0)
MTSV<-left_join(MTSV,addmich)

fruiting<-data.frame(name=MTSV$name,fruiting=MTSV$fruiting)
HF.data<-left_join(HF.data,fruiting)
HF.data<-left_join(HF.data,addmich)

library(xtable)
?xtable()
floweringcors<-data.frame(data=c("MTSV","HF"),fruit.development=c(cor(MTSV$flo_time,MTSV$fruiting,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$fruiting,use="pairwise.complete.obs")),
 seed.mass=c(cor(MTSV$flo_time,MTSV$seed_mass,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$seed_mass,use="pairwise.complete.obs")),
 xylem_anatomy=c(cor(MTSV$flo_time,MTSV$xylem_anatomy,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$xylem_anatomy,use="pairwise.complete.obs")),
  cold.tol=c(cor(MTSV$flo_time,MTSV$min_temp,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$min_temp,use="pairwise.complete.obs")),
             growing.season=c(cor(MTSV$flo_time,MTSV$frost_free,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$frost_free,use="pairwise.complete.obs")),
           bud.type=c(cor(MTSV$flo_time,MTSV$flo_loc_num,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$flo_loc_num,use="pairwise.complete.obs")))

xtable(floweringcors,digits=3)


MTSV.drought<-left_join(MTSV,freezish)



minnyP<-data.frame(predictor=c("min.precip"),ψ50=cor(MTSV.drought$min._precip,MTSV.drought$ψ50,use="pairwise.complete.obs"),
           ψ50.safety.margin=cor(MTSV.drought$min._precip,MTSV.drought$ψψ50.safety.margin,use="pairwise.complete.obs"),
           Aridity.index=cor(MTSV.drought$min._precip,MTSV.drought$AI..aridity.index.,use="pairwise.complete.obs"),
           cold.tol=cor(MTSV.drought$min._precip,MTSV.drought$min_temp,use="pairwise.complete.obs"),
           xylem.anatomy=cor(MTSV.drought$min._precip,MTSV.drought$xylem_anatomy,use="pairwise.complete.obs"))

xtable(minnyP,digits=3)

data.frame(predictor=c("min_P"),ψ50=cor(MTSV.drought$min._precip,MTSV.drought$ψ50,use="pairwise.complete.obs")
,           xylem.anatomy=cor(MTSV.drought$xylem_anatomy,MTSV.drought$ψ50,use="pairwise.complete.obs"),
           flotime=cor(MTSV.drought$flo_time,MTSV.drought$ψ50,use="pairwise.complete.obs") )




data.frame(predictor=c("P50"),cold.tol=cor(MTSV.drought$min_temp,MTSV.drought$ψ50,use="pairwise.complete.obs"),
          xylem.anatomy=cor(MTSV.drought$xylem_anatomy,MTSV.drought$ψ50,use="pairwise.complete.obs"),
          flotime=cor(MTSV.drought$flo_time,MTSV.drought$ψ50,use="pairwise.complete.obs") )
           

cor(MTSV.drought$seed_mass,MTSV.drought$flo_time,use="pairwise.complete.obs")
MTSV.drought$fruit_size<-as.numeric(MTSV.drought$fruit_size)

cor(MTSV.drought$flo_time,MTSV.drought$fruit_size,use="pairwise.complete.obs")
cor(MTSV.drought$flo_time,MTSV.drought$min_temp,use="pairwise.complete.obs")
cor(MTSV.drought$flo_time,MTSV.drought$frost_free,use="pairwise.complete.obs")
cor(MTSV.drought$av_fruit_time,MTSV.drought$frost_free,use="pairwise.complete.obs")
cor(MTSV.drought$av_fruit_time,MTSV.drought$flo_time,use="pairwise.complete.obs")
