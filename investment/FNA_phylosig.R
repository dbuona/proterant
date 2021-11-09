### Started by Dan Feb 9 2021
### Part 1:### modeling FLS for FNA with measurement models in brms
#FLS ~ pdsi + flower traits + fruit traits
###Part 2: Focusing on prunocerasus with stan model
#based on 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
options(mc.cores = parallel::detectCores())


library(dplyr)
library(tidybayes)
library(bayesplot)
library("caper")
library(phytools)
library(geiger)


#-------------------------#
#------Part 1--------------#
#--------------------------#
setwd("~/Documents/git/proterant/investment")

FNA<-read.csv("Data/cherry_data.csv") ## measurement and FLS data from FNA
##clean specicies name

FNA$species<-ifelse(FNA$species=="hortunlana","hortulana",FNA$species)
FNA$species<-ifelse(FNA$species=="gladulosa","glandulosa",FNA$species)
FNA$species<-ifelse(FNA$species=="speciosa","speciosa",FNA$species)
FNA$species<-ifelse(FNA$species=="fasiculata","fasciculata",FNA$species)

FNA$FLSnum<-NA
FNA$FLSnum[FNA$FLS=="before"]<-1
FNA$FLSnum[FNA$FLS=="before/with"]<-2
FNA$FLSnum[FNA$FLS=="with"]<-3
FNA$FLSnum[FNA$FLS=="after"]<-4
FNA$FLSnum<-as.integer(FNA$FLSnum)

##read in tree, root it and give it branch lengths
mytree<-read.tree(file = "cherrytree.tre")
plotTree(mytree)
mytree<-root(mytree,outgroup = "Physocarpus")
mytree2<-compute.brlen(mytree, method = "Grafen")
plotTree(mytree2)
is.ultrametric(mytree2)

####now need to subset data
names.intree<-mytree2$tip.label 

FNA.small<-dplyr::filter(FNA,species %in% names.intree)

setdiff(FNA$species,mytree2$tip.label)
setdiff(mytree2$tip.label,FNA$species)

setdiff(FNA.small$species,mytree2$tip.label)
setdiff(mytree2$tip.label,FNA.small$species)

namelist<-unique(FNA.small$species)

to.prune<-which(!names.intree%in%namelist) #prun the tree
pruned.by<-drop.tip(mytree2,to.prune)
plotTree(pruned.by)

setdiff(FNA.small$species,pruned.by$tip.label)
setdiff(pruned.by$tip.label,FNA.small$species)

mytree.names<-pruned.by$tip.label 
#pruned.by$node.label<-NULL




length(FNA.small$species)
length(pruned.by$tip.label)

FNA.small2<-dplyr::select(FNA.small,species,FLSnum)
final.df<-FNA.small2[match(mytree.names, FNA.small2$species),]
namelist2<-final.df$species
namelist2==mytree.names
final.df$species== mytree.names





phylosig(pruned.by,final.df$FLSnum,method="lambda",nsim = 100, test=TRUE)

library(ggtree)
mytree2$tip.labels<-final.df$species
full_join(mytree2,final.df)
final.df$FLSfact<-as.factor(final.df$FLSnum)


jpeg("Plots/phylosig1.jpeg",height=6,width=6,units='in',res=300)
p<-ggtree(pruned.by)
p %<+% final.df+geom_tiplab(hjust=-.2,align=TRUE,fontface="italic")+geom_tippoint(aes(color=FLSfact),size=5)+ xlim(0, 2)+geom_cladelabel(node=4,offset=-.7, label=" lambda= 0.34")+scale_color_viridis_d(option="plasma",name="Flower-leaf sequence",labels=c("before","before/with","with","after"))
dev.off()
?scale_color_viridis_d()
fitLambda<-function(tree,x,model="ER"){
  lik<-function(lambda,tree,x,model)
    logLik(ace(x,rescale(tree,model="lambda",lambda),
               type="discrete",model=model))
  obj<-optimize(lik,c(0,1),tree=tree,x=x,model=model,maximum=TRUE)
  fit<-ace(x,rescale(tree,model="lambda",lambda=obj$maximum),
           type="discrete",model=model)
  I<-fit$index.matrix
  fitted.Q=matrix(fit$rates[I],dim(I)[1],dim(I)[2],
                  dimnames=list(dimnames(fit$lik.anc)[[2]],
                                dimnames(fit$lik.anc)[[2]]))
  diag(fitted.Q)<--rowSums(fitted.Q,na.rm=TRUE)
  list(Q=fitted.Q,lambda=obj$maximum,logLik=logLik(fit))
}
fitLambda(pruned.by,final.df$FLSnum)

pruned.by
### Q1 DO we bother with phylogeny? Only 7 species, and as you can see, weak signal

### phylogeny model
write.tree(pruned.by,"pruned_4_modeling.tre")


