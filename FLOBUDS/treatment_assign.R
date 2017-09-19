##This code uses fake data to simulate the treatment assignment process in a block.
setwd("~/Documents/git/proterant/FLOBUDS")
R<-read.csv("randomize.csv",header=TRUE)
library(randomizr)
N<-nrow(R)
Z <- complete_ra(N = N)
table(Z)


Z <- complete_ra(N = N, num_arms = 8)
table(Z)


### I think its like this but with species
Z <- block_ra(block_var = R$species, num_arms = 8)
R$assignment<-Z

Ac<-filter(R,genus=="Acer")
table(Ac$assignment)

Bet<-filter(R,genus=="Betula")
table(Bet$assignment)

Cor<-filter(R,genus=="Corylus")
table(Cor$assignment)

### asign individual to beaker
L<-filter(R, assignment=="T1")

Z <- complete_ra(N = nrow(L), num_arms = 2)


