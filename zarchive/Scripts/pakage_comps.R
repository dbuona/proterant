rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library("MCMCglmm")
library(phylolm)
library(brms)

n <- 100
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

# Generate random data and standardize to have mean 0 and variance 1
X1 <- rTraitCont(phy, model = "BM", sigma = 1)
X1 <- (X1 - mean(X1))/var(X1)
X2 <- rTraitCont(phy, model = "BM", sigma = 1)
X2 <- (X2 - mean(X1))/var(X2)



# Simulate binary Y
sim.dat <- data.frame(Y=array(0, dim=n), X1=X1,X2=X2, row.names=phy$tip.label)

sim.dat$Y <- binaryPGLMM.sim(Y ~ X1+X2, phy=phy, data=sim.dat, s2=0,
                             B=matrix(c(1,4,0.5),nrow=3,ncol=1), nrep=1)$Y

# Fit model
summary(glm(Y ~ X1+X2,data = sim.dat,family="binomial"))
library("logistf")
logistf(Y ~ X1+X2, data=sim.dat)

binaryPGLMM(Y ~ X1+X2, phy=phy, data=sim.dat)
summary(phyloglm(Y ~ X1+X2, phy=phy, data=sim.dat,method = "logistic_MPLE", btol = 200, log.alpha.bound = 6))





sim.dat$species <- phy$tip.label
inv.phylo <- MCMCglmm::inverseA(phy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


get_prior(Y~X1+(1|species), data =sim.dat, family = bernoulli(link="logit"))

a<-brm(Y~X1+(1|species), data =sim.dat, family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000,thin=1,  prior = c(
  prior(normal(1,3), "b"),
  prior(student_t(3, 0, 10), "Intercept"),
  prior(student_t(3, 0, 10), "sd")))
summary(a)### brms inflates predic
pp_check(a)

hist(rnorm(1000,1,3))
