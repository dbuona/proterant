setwd("~/Desktop/")

library(ape)

a<-read.tree("..//Downloads/sphingidae.randomly.resolved.tre")

plot(a)
a

write.tree(a,"..//Documents/git/proterant/investment/Input/plum.tre")
