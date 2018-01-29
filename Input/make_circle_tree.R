####This makes a phylogenetic tree with hysteranthy as a trait. 1/10/18

#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")


library(ggtree)


again<-read.csv("mich_data_full.csv")

tr <-mich.tree

par(mar=c(.2,.2,.2,.2))
p<-ggtree(tr,layout="circular")



dd <- data.frame(taxa  = again$name, hysteranthy = as.character(again$pro2 ))
row.names(dd) <- NULL
print(dd)

p <- p %<+% dd + geom_tiplab2(size=2)+geom_tippoint(aes(color=hysteranthy,shape=hysteranthy))
p+theme(legend.position="right")

