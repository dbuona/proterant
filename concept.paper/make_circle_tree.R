####This makes a phylogenetic tree with hysteranthy as a trait. 1/10/18

#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")


library(ggtree)
setwd("~/Documents/git/proterant/input")


mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
mich.data$synth<-mich.data$pro2+mich.data$pro3
mich.data$synth2[which(mich.data$synth==0)] <- "always seranthous"
mich.data$synth2[which(mich.data$synth==1)] <- "transitional"
mich.data$synth2[which(mich.data$synth==2)] <- "always hysteranthous"

tr <-mich.tree

par(mar=c(.1,.1,.1,.1))
p<-ggtree(tr,layout="circular")

dd <- data.frame(taxa  = mich.data$name, hysteranthy = mich.data$synth2 )
#row.names(dd) <- NULL
print(dd)

p <- p %<+% dd + geom_tiplab2(size=1.8)+geom_tippoint(aes(color=hysteranthy,shape=hysteranthy))
p+theme(legend.position="right")


silv.data$synth<-silv.data$pro2+silv.data$pro3
silv.data$synth2[which(silv.data$synth==0)] <- "always seranthous"
silv.data$synth2[which(silv.data$synth==1)] <- "transitional"
silv.data$synth2[which(silv.data$synth==2)] <- "always hysteranthous"

#### make a tree for svics

tr1 <-silv.tree
q<-ggtree(tr1,layout="circular")
dd1 <- data.frame(taxa  = silv.data$name, hysteranthy = silv.data$synth2 )

q <- q %<+% dd1 + geom_tiplab2(size=1.8)+geom_tippoint(aes(color=hysteranthy,shape=hysteranthy))
q+theme(legend.position="right")
