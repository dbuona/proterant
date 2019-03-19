rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")

library(MCMCglmm)
#if you dont want to run the model: 
load("RData/zarchival/hystmodels.RData")

#1) need to get a VCV from the phylgeny

#michVCV<-vcv(mich.tree.droughtprune)
phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
data_simple <- read.table(
  "https://paul-buerkner.github.io/data/data_simple.txt", 
  header = TRUE
)
head(data_simple)


##construct covarience matrix:
mich.data<-rownames_to_column(mich.data, "name")
inv.phylo <- MCMCglmm::inverseA(mich.tree.droughtprune, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

A
###bayesian and continuous-- main model###############
modelcheck <- brm(pro2~ pol_cent+flo_cent+precip_cent+ (1|name), data = mich.data, 
                 family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000)
                 prior = c(prior(normal(0, 5), "b"),
                          prior(normal(0, 5), "Intercept"),
                           prior(student_t(3, 0, 5), "sd"))) 

summary(modelcheck)

colnames(mich.data)
MacMac <- MCMCglmm(pro2~ pol_cent+flo_cent+precip_cent,random=~name, family ="ordinal",data=mich.data, ginverse = list(name=A))
summary(MacMac)
?MCMCglmm()
summary(z.funct.drought.noint)

Lmat <- matrix(rep(1), nrow = nrow(michVCV), ncol = ncol(michVCV))
diag(Lmat) <- 0

datalist<- with(mich.data, 
                list(y=pro2,
                     pol=pol_cent,
                     flotime=flo_cent,
                     minP=precip_cent,
                     V=michVCV,
                     Lmat=Lmat,
                     N = nrow(mich.data)
                ))

z.funct.stan.phylo<- stan('..//stan/pgls_hyst2.stan', data = datalist,
             iter = 3000, warmup=2000) 
