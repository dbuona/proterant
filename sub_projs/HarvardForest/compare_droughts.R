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

data.frame(data=c("MTSV","HF"),fruit.development=c(cor(MTSV$flo_time,MTSV$fruiting,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$fruiting,use="pairwise.complete.obs")),
 seed.mass=c(cor(MTSV$flo_time,MTSV$seed_mass,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$seed_mass,use="pairwise.complete.obs")),
 xylem_anatomy=c(cor(MTSV$flo_time,MTSV$xylem_anatomy,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$xylem_anatomy,use="pairwise.complete.obs")),
  cold.tol=c(cor(MTSV$flo_time,MTSV$min_temp,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$min_temp,use="pairwise.complete.obs")),
             growing.season=c(cor(MTSV$flo_time,MTSV$frost_free,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$frost_free,use="pairwise.complete.obs")),
           bud.type=c(cor(MTSV$flo_time,MTSV$flo_loc_num,use="pairwise.complete.obs"),cor(HF.data$fopn.jd,HF.data$flo_loc_num,use="pairwise.complete.obs")))



MTSV.drought<-left_join(MTSV,freezish)

data.frame(predictor=c("min.precip"),ψ50=cor(MTSV.drought$min._precip,MTSV.drought$ψ50,use="pairwise.complete.obs"),
           ψ50.safety.margin=cor(MTSV.drought$min._precip,MTSV.drought$ψψ50.safety.margin,use="pairwise.complete.obs"),
           Aridity.index=cor(MTSV.drought$min._precip,MTSV.drought$AI..aridity.index.,use="pairwise.complete.obs"),
           cold.tol=cor(MTSV.drought$min._precip,MTSV.drought$min_temp,use="pairwise.complete.obs"),
           xylem.anatomy=cor(MTSV.drought$min._precip,MTSV.drought$xylem_anatomy,use="pairwise.complete.obs"))


data.frame(predictor=c("P50"),cold.tol=cor(MTSV.drought$min_temp,MTSV.drought$ψ50,use="pairwise.complete.obs"),
          xylem.anatomy=cor(MTSV.drought$xylem_anatomy,MTSV.drought$ψ50,use="pairwise.complete.obs"),
          flotime=cor(MTSV.drought$flo_time,MTSV.drought$ψ50,use="pairwise.complete.obs") )
           
HFtree<-read.tree("HarvardForest/HFtree4modeling.tre" )

####mean only model
HF.tree<-drop.tip(HF.tree,tip = "Viburnum_lantanoides")
HF.tree$tip.label

HF.means<-HF.data %>% group_by(name)%>% summarise(meanFLS.func=mean(funct.fls,na.rm=TRUE),meanFLS.phys=mean(phys.fls,na.rm=TRUE),meanFLS.inter=mean(inter.fls,na.rm=TRUE),meanflo=mean(flo_cent,na.rm=TRUE),meandrought=mean(precip_cent),meanpol=mean(pol_cent))
HF.means$name

HF.means$hyster.fuct<-ifelse(HF.means$meanFLS.func>1,1,0)
HF.means$hyster.phys<-ifelse(HF.means$meanFLS.phys>1,1,0)
HF.means$hyster.inter<-ifelse(HF.means$meanFLS.inter>1,1,0)
HF.means<-  HF.means %>% remove_rownames %>% column_to_rownames(var="name")

z.phys.HFmean<-phylolm(meanFLS.phys~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean)

z.phys.HFbin<-phylolm(hyster.phys~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

z.funct.HFmean<-phylolm(meanFLS.func~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.funct.HFmean)

z.funct.HFbin<-phylolm(hyster.fuct~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)


z.inter.HFmean<-phylolm(meanFLS.inter~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree,model="BM", boot=599,full.matrix = TRUE)
summary(z.inter.HFmean)

z.inter.HFbin<-phylolm(hyster.inter~meanflo+meandrought+meanpol+meanflo:meandrought+meanflo:meanpol+meanpol:meandrought,HF.means,HF.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=599,full.matrix = TRUE)

summary(z.funct.HFbin)
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
}
extract_CIs<-function(x){
  dplyr::filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
}

HF.funct.quan<-full_join(extract_coefs(z.funct.HFmean),extract_CIs(z.funct.HFmean),by="trait")

colnames(HF.funct.quan)<-c("trait","estimate","low","high")
HF.funct.quan$class<-"functional"

HF.phys.quan<-full_join(extract_coefs(z.phys.HFmean),extract_CIs(z.phys.HFmean),by="trait")
colnames(HF.phys.quan)<-c("trait","estimate","low","high")
HF.phys.quan$class<-"physiological"

quanty<-rbind(HF.phys.quan,HF.funct.quan)
quanty$model<-"quantitative" 

HF.funct.bin<-full_join(extract_coefs(z.funct.HFbin),extract_CIs(z.funct.HFbin),by="trait")

colnames(HF.funct.bin)<-c("trait","estimate","low","high")
HF.funct.bin$class<-"functional"

HF.phys.bin<-full_join(extract_coefs(z.phys.HFbin),extract_CIs(z.phys.HFbin),by="trait")
colnames(HF.phys.bin)<-c("trait","estimate","low","high")
HF.phys.bin$class<-"physiological"

binny<-rbind(HF.phys.bin,HF.funct.bin)
binny$model<-"categorical" 

comps<-rbind(binny,quanty)
comps<-dplyr::filter(comps,trait!="(Intercept)")
comps<-dplyr::filter(comps,trait!="sigma2")
unique(comps$trait)
comps$trait[which(comps$trait=="meanpol")]<-"pollination syndrome"
comps$trait[which(comps$trait=="meanflo")]<- "early flowering"
comps$trait[which(comps$trait=="meandrought")]  <- "water dynamics"
comps$trait[which(comps$trait=="meandrought:meanpol")]<- "pollination:water dynamics"
comps$trait[which(comps$trait=="meanflo:meanpol")]<-"pollination:early flowering"
comps$trait[which(comps$trait=="meanflo:meandrought")]<-"early flowering:water dynamics"
pd=position_dodgev(height=0.4)

physcomps<-filter(comps,class=="physiological")
functscomps<-filter(comps,class=="functional")

physcomps$class<-"fbb-lbb"
physcomps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class,color=model),position=pd,size=3)+
  geom_errorbarh(aes(xmin=low,xmax=high,linetype=class,color=model),position=pd,height=0,size=.5)+
  scale_linetype_manual(values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  scale_color_manual(values=c("orchid4","springgreen4"))+guides(size = "legend", linetype= "none")+facet_wrap(~model,scales = "free_x")

functscomps$class<-"fopn-l75"
functscomps %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("early flowering:water dynamics","pollination:early flowering","pollination:water dynamics","early flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=class,color=model),position=pd,size=3)+
  geom_errorbarh(aes(xmin=low,xmax=high,linetype=class,color=model),position=pd,height=0,size=.5)+
  scale_linetype_manual(values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  scale_color_manual(values=c("orchid4","springgreen4"))+guides(size = "legend", linetype= "none")+facet_wrap(~model,scales = "free_x")

