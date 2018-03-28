rm(list=ls())
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(randomForest)
library(gbm)

library(pez)

setwd("C:/Users/Ignacio/Documents/MEGA/Work_Montreal_postdoc/sWEEP/Joanne Project/TTOL data")
load(".RData")


ttol<-read.newick("2.TTOL_all_smoothed.nwk")
ttol2<-read.newick("1.TTOL_all_unsmoothed.nwk")


names.inTTOL<-ttol$tip.label
names.inTTOL2<-ttol2$tip.label

genus.inTTOL<-sapply(strsplit(names.inTTOL, "_"), "[", 1)
genus.inTTOL2<-sapply(strsplit(names.inTTOL2, "_"), "[", 1)

unique.genus.TTOL<-unique(genus.inTTOL)
unique.genus.TTOL2<-unique(genus.inTTOL2)

unique.genus.SWEEP<-unique(sWEEP.data$Genus)

represented.genus<-which(unique.genus.SWEEP%in%unique.genus.TTOL)
represented.genus2<-which(unique.genus.SWEEP%in%unique.genus.TTOL2)

genus.to.prune<-which(!genus.inTTOL%in%unique.genus.SWEEP)
genus.to.prune2<-which(!genus.inTTOL2%in%unique.genus.SWEEP)

ttol.pruned.by.genus<-drop.tip(ttol,genus.to.prune)
ttol2.pruned.by.genus<-drop.tip(ttol2,genus.to.prune2)

names.inTTOL.pruned<-ttol.pruned.by.genus$tip.label
names.inTTOL.pruned2<-ttol2.pruned.by.genus$tip.label

genus.inTTOL.pruned<-sapply(strsplit(names.inTTOL.pruned, "_"), "[", 1)
genus.inTTOL.pruned2<-unique(sapply(strsplit(names.inTTOL.pruned2, "_"), "[", 1))


sweep.compound.names<-paste(sWEEP.data$Genus,sWEEP.data$Species,sep="_")
sweep.compound.names2<-paste(sWEEP.data$Genus,sWEEP.data$Species,sep="_")
sps.to.prune<-which(!ttol.pruned.by.genus$tip.label%in%sweep.compound.names)
sps.to.prune2<-which(!ttol2.pruned.by.genus$tip.label%in%sweep.compound.names2)


ttol.sweep.phy<-drop.tip(ttol.pruned.by.genus,sps.to.prune)
ttol.sweep.phy2<-drop.tip(ttol2.pruned.by.genus,sps.to.prune2)



########################################################################################################
####################################################
##
##        final phylogenies
##
####################################################
########################################################################################################



## loading data
sWEEP.final<-read.csv(file="SWEEP_DATA_20_07.csv", header=T,sep=",")

head(sWEEP.final)
colnames(sWEEP.final)



############## Defining data
sWEEP.data.cor<-read.table("SWEEP_DATA_MS_1_ages2.txt",as.is=T,header=T,sep="\t")
names(sWEEP.data.cor)
sWEEP.data<-sWEEP.data.cor
sWEEP.data$thermy<-ifelse(sWEEP.data$Class %in% c("Mammalia", "Aves"), "endo", "ecto")
sWEEP.data$is.aquatic<-ifelse(sWEEP.data$Realm %in% c("Marine", "Freshwater"), "Aquatic", 
                              ifelse(sWEEP.data$Realm %in% "Intertidal",NA,"Terrestrial"))
sWEEP.data$Palaeotemp<-ifelse(sWEEP.data$order.temp.origin ==5, 3, 
                              ifelse(sWEEP.data$order.temp.origin==4,2,1))
sWEEP.data$warm.cold<-ifelse(sWEEP.data$Palaeotemp ==1, 1,2) 

sWEEP.data$Plant<-ifelse(sWEEP.data$Phylum %in%c("Phaeophyceae",
                                                 "Streptophyta","Ascomycota","Basidiomycota"),
                         "Plants","No.plants")

ectolandwarmo<-subset(sWEEP.data, sWEEP.data$thermy=="ecto"& sWEEP.data$Realm=="Terrestrial" & sWEEP.data$order.temp.origin==5)

sWEEP.final<-sWEEP.data

## species level phylogeny
sWEEP.final$species.phylo<-paste(sWEEP.final$Genus,sWEEP.final$Species,sep="_")

to.rem.ttolsps<-which(!ttol$tip.label%in%sWEEP.final$species.phylo)
phy.sweep.sps<-drop.tip(ttol,to.rem.ttolsps)

to.rem.ttolsps<-which(!ttol2$tip.label%in%sWEEP.final$species.phylo)
phy.sweep.sps<-drop.tip(ttol2,to.rem.ttolsps)


#### species level phylo-signal
parslist<-list(c("ecto","No.plants"),
               c("endo","No.plants"),
               c("ecto","Plants"),
               c("ecto",1,"No.plants"),
               c("ecto",2,"No.plants"),
               c("endo",1,"No.plants"),
               c("endo",2,"No.plants"))

storing.phylosig<-as.data.frame(array(NA,dim=c(8,6)))
for(i in 4:7){
  print(i)
  phy.t<-phy.sweep.sps
  
  sWEEP.final<-subset(sWEEP.data,sWEEP.data$thermy==parslist[[i]][1]
                      & sWEEP.data$warm.cold==parslist[[i]][2]
                      & sWEEP.data$Plant==parslist[[i]][3])
  sWEEP.final$species.phylo<-paste(sWEEP.final$Genus,sWEEP.final$Species,sep="_")
  
  
  x.up<-sWEEP.final$tmax
  x.low<-sWEEP.final$tmin
  
  sWEEP.final<-sWEEP.data
  names(x.up)<-sWEEP.final$species.phylo
  names(x.low)<-sWEEP.final$species.phylo
  
  # remove nas
  #x.up<-x.up[which(!is.na(x.up))]
  #x.low<-x.low[which(!is.na(x.low))]
  
  ## create comparative data
  
  comp.sps.tmax<-comparative.data(phy.t,data.frame(x.up,namestmax=names(x.up)),names.col="namestmax")
  comp.sps.tmin<-comparative.data(phy.t,data.frame(x.low,namestlow=names(x.low)),names.col="namestlow")
  
  
  ## remove non in data
  phy.low<-comp.sps.tmin$phy
  phy.up<-comp.sps.tmax$phy
  
  #rownames(comp.sps.tmax$data)==comp.sps.tmax$phy$tip.label
  vec.up<-as.vector(comp.sps.tmax$data[,1])
  vec.low<-as.vector(comp.sps.tmin$data[,1])
  names(vec.up)<-rownames(comp.sps.tmax$data)
  names(vec.low)<-rownames(comp.sps.tmin$data)
  
  #up.delta<-fitContinuous(comp.sps.tmax$phy,x.up[match(comp.sps.tmax$phy$tip.label,names(x.up))],model="delta",bounds=list(delta=c(0.1,100)))
  #low.delta<-fitContinuous(comp.sps.tmin$phy,x.low[match(comp.sps.tmin$phy$tip.label,names(x.low))],model="delta",bounds=list(delta=c(0.1,100)))
  phy.upOU<-rescaleTree(comp.sps.tmax$phy,1)
  phy.lowOU<-rescaleTree(comp.sps.tmin$phy,1)
  up.OU<-fitContinuous(phy.upOU,vec.up,model="OU",bounds=list(alpha=c(0.0005,20)))
  low.OU<-fitContinuous(phy.lowOU,vec.low,model="OU",bounds=list(alpha=c(0.0005,20)))
  up.BM<-fitContinuous(comp.sps.tmax$phy,vec.up,model="BM")
  low.BM<-fitContinuous(comp.sps.tmin$phy,vec.low,model="BM")
  up.lambda<-fitContinuous(comp.sps.tmax$phy,vec.up,model="white")
  low.lambda<-fitContinuous(comp.sps.tmin$phy,vec.low,model="white")
  ?fitContinuous
  
  storing.phylosig[impar[i],1]=length(comp.sps.tmax$phy$tip.label)
  storing.phylosig[impar[i],2]=up.BM$opt$sigsq
  storing.phylosig[impar[i],3]=log(up.OU$opt$alpha,10)*(-1)
  storing.phylosig[impar[i],4]=up.OU$opt$lnL
  storing.phylosig[impar[i],5]=up.BM$opt$lnL
  storing.phylosig[impar[i],6]=up.lambda$opt$lnL
  storing.phylosig[par[i],1]=length(comp.sps.tmin$phy$tip.label)
  storing.phylosig[par[i],2]=low.BM$opt$sigsq
  storing.phylosig[par[i],3]=log(low.OU$opt$alpha,10)*(-1)
  storing.phylosig[par[i],4]=low.OU$opt$lnL
  storing.phylosig[par[i],5]=low.BM$opt$lnL
  storing.phylosig[par[i],6]=low.lambda$opt$lnL
  
}

write.table(storing.phylosig,file="phylosig.endo_ecto.cold.warm.unsmoothed.txt",sep="\t",col.names = T,row.names=T)



######### MODE OF EVOLUTION - WHAT TRAITGRAMS SHOW
library(pez)
library(picante)
library(phytools)
library(caper)
#phylo and traits
phy.t<-phy.sweep.sps

sWEEP.final.t<-subset(sWEEP.data,sWEEP.data$thermy==parslist[[1]][1]
                    #& sWEEP.data$warm.cold==parslist[[]][2]
                    & sWEEP.data$Plant==parslist[[3]][2])


sWEEP.final.t$species.phylo<-paste(sWEEP.final.t$Genus,sWEEP.final.t$Species,sep="_")
x.up<-sWEEP.final.t$tmax
x.low<-sWEEP.final.t$tmin
names(x.up)<-sWEEP.final.t$species.phylo
names(x.low)<-sWEEP.final.t$species.phylo


## create comparative data
comp.sps.tmax<-comparative.data(phy.t,data.frame(x.up,namestmax=names(x.up)),names.col="namestmax")
comp.sps.tmin<-comparative.data(phy.t,data.frame(x.low,namestlow=names(x.low)),names.col="namestlow")

fff<-traitgram.cc(comp.sps.tmax,comp.sps.tmax$data$x.up)

vec<-comp.sps.tmax$data$x.up
names(vec)<-comp.sps.tmax$phy$tip.label
phenogram(comp.sps.tmax$phy,vec)
phenogram(comp.sps.tmax$phy,vec, spread.labels = TRUE, spread.cost = c(1, 0))
fancyTree(comp.sps.tmax$phy, type = "phenogram95", x = vec, spread.cost = c(1, 0))





############# RANDOM FORESTS


############## Defining data
sWEEP.data.cor<-read.table("SWEEP_DATA_MS_1_ages2.txt",as.is=T,header=T,sep="\t")
names(sWEEP.data.cor)
dim(sWEEP.data)

ectolandwarmo<-subset(sWEEP.data, sWEEP.data$thermy=="ecto"& sWEEP.data$Realm=="Terrestrial" & sWEEP.data$order.temp.origin==5)

sWEEP.final<-sWEEP.data

store.random.forest<-array(NA,dim=c(6,5))
parslist<-list(c("ecto","No.plants"),
               c("endo","No.plants"),
               c("ecto","Plants"))

for(i in 1:3){
  print(i)
  subs.data<-subset(sWEEP.data,sWEEP.data$thermy==parslist[[i]][1]&sWEEP.data$Plant==parslist[[i]][2])
  subs.data$order.age<-as.numeric(scale(subs.data$order.age))
  subs.data$order.temp.origin<-as.numeric(scale(subs.data$order.temp.origin))
  subs.data$lat_max<-as.numeric(scale(as.numeric(subs.data$lat_max)))
  subs.data$lat_min<-as.numeric(scale(as.numeric(subs.data$lat_min)))
  
  mod.tmax1<-randomForest(tmax~order.age+order.temp.origin+lat_max,data=subs.data,
                          ntree=500,na.action=na.omit,importance=T,nPerm=50,mtry=2,corr.bias=T,localImp=T)
  mod.tmin1<-randomForest(tmin~order.age+order.temp.origin+lat_min,data=subs.data,
                          ntree=500,na.action=na.omit,importance=T,nPerm=50,mtry=2,corr.bias=T,localImp=T)
  #mod.tmax1<-randomForest(tmax~order.age+order.temp.origin,data=subs.data,
   #                       ntree=500,na.action=na.omit,importance=T,nPerm=50,mtry=2,corr.bias=T,localImp=T)
  #mod.tmin1<-randomForest(tmin~order.age+order.temp.origin,data=subs.data,
   #                       ntree=500,na.action=na.omit,importance=T,nPerm=50,mtry=2,corr.bias=T,localImp=T)
  
 
    store.random.forest[impar[i],1]<-mod.tmax1$rsq[500]
    store.random.forest[par[i],1]<-mod.tmin1$rsq[500]
    
    store.random.forest[impar[i],2]<-mod.tmax1$mse[500]
    store.random.forest[par[i],2]<-mod.tmin1$mse[500]
    
    store.random.forest[impar[i],3:5]<-mod.tmax1$importance[,1]
    store.random.forest[par[i],3:5]<-mod.tmin1$importance[,1]

      
      }

write.table(store.random.forest,file="random.forest.full2.txt",sep="\t",col.names = T,row.names=T)
#write.table(store.random.forest,file="random.forest.no.lats.txt",sep="\t",col.names = T,row.names=T)

## correlations
list.cors<-list()
for(i in 1:3){
subs.data<-subset(sWEEP.data,sWEEP.data$thermy==parslist[[i]][1]&sWEEP.data$Plant==parslist[[i]][2])
subs.data$lat_max<-as.numeric(subs.data$lat_max)
subs.data$lat_min<-as.numeric(subs.data$lat_min)
list.cors[[impar[i]]]<-cor(cbind(subs.data$tmax,subs.data$order.age,subs.data$order.temp.origin,subs.data$lat_max),use="pairwise.complete.obs")
list.cors[[par[i]]]<-cor(cbind(subs.data$tmin,subs.data$order.age,subs.data$order.temp.origin,subs.data$lat_max),use="pairwise.complete.obs")
}




## plotting CTAs
# Regression Tree Example
#library(rpart)
#library(rpart.plot)
library(party)

## prepare data and fit trees
i<-2
sWEEP.data<-sWEEP.data.cor
subs.data<-subset(sWEEP.data,sWEEP.data$thermy==parslist[[i]][1]&sWEEP.data$Plant==parslist[[i]][2])
subs.data<-subset(subs.data,!is.na(subs.data$tmax))
subs.data$lat_max<-as.numeric(subs.data$lat_max)
fit.tmax <- ctree(tmax~order.age+order.temp.origin+lat_max,data=subs.data
            ,controls=ctree_control(maxdepth =3,mincriterion = 0.7))
plot(fit.tmax)


subs.data<-subset(sWEEP.data,sWEEP.data$thermy==parslist[[i]][1]&sWEEP.data$Plant==parslist[[i]][2])
subs.data<-subset(subs.data,!is.na(subs.data$tmin))
subs.data$lat_min<-as.numeric(subs.data$lat_min)
fit.tmin <- ctree(tmin~order.age+order.temp.origin+lat_min,data=subs.data
            ,controls=ctree_control(maxdepth =3,mincriterion = 0.55))

plot(fit.tmin)



## alternative method
tree <- rpart(tmax~order.age+order.temp.origin+lat_max,data=subs.data, method="anova"
              ,control=rpart.control(maxdepth =3))
#,control=rpart.control(cp =0.05))
tree <- rpart(tmax~order.age+order.temp.origin,data=subs.data, method="anova"
              ,control=rpart.control(maxdepth =4))
#,control=rpart.control(cp =0.05))

