###This is Dan's main thesis file.   

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library(ape)
library(phytools)
library(geiger)
library(gbm)
library(pez)
library(dplyr)
library(tidyr)
library(caper)
library(picante)
library(tidyverse)
library(boot)
library(phylolm)
#https://academic-oup-com.ezp-prod1.hul.harvard.edu/sysbio/article-lookup/doi/10.1093/sysbio/syp074 Garland and Ives 2010

###read in tree from Zanne et al
treee<-read.tree("Vascular_Plants_rooted.dated.tre")
is.ultrametric(treee)### is not ultrametric
anthy<-read.csv("michigantrees_sequence.csv", header = TRUE)
anthy<-filter(anthy, !is.na(av_fruit_time))
source("source/prune_tree.R")
is.ultrametric(pruned.by.anthy)
#write.tree(pruned.by.anthy, "michigan.phy")
########################################################################
#phlogenetic signal####################################################
##WOrked! but only when the data dpoesnt match

###try making it a comaprative data object for phylo.D
d<-comparative.data(pruned.by.anthy,final.df,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)

plot(pruned.by.anthy)
PhyloD <- phylo.d(d, binvar=pro)##D=0.3156129
summary(pruned.by.anthy)


names<-final.df$name
tips<-pruned.by.anthy$tip.label
names==tips

###now is pollination still important for within early ones
#data
earl<-filter(final.df, fruit_bin==0)
## tree
namelist2<-unique(earl$name)
names.intree2<-pruned.by.anthy$tip.label
##Prune the tree
to.prune<-which(!names.intree2%in%namelist2)
pruned.earl<-drop.tip(pruned.by.anthy,to.prune)

##MODELS#########################################################################

#make $name row names
final.df<-  final.df %>% remove_rownames %>% column_to_rownames(var="name")


####full model with everything bianary###
full.mod<-glm(pro~pol+class2+shade_bin+fruit_bin+flo_type,family = binomial(link="logit"),data=final.df)
summary(full.mod)

####full  phylogentically corrected########### this model seems sensative when to the random additions? sometimes shade is significant 
full.modA<-phyloglm(pro~pol+class2+shade_bin+fruit_bin+flo_type,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                start.beta=NULL, start.alpha=NULL,
                 boot=10,full.matrix = TRUE)
summary(full.modA)


###model with fruit time and height as continuous, addind flower time
full.modB<-phyloglm(pro~pol+heigh_height+shade_bin+av_fruit_time+flo_time+flo_type,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 100, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=10,full.matrix = TRUE)
summary(full.modB) ### signifcance and direction does not change with coninuous, height does become signifiant marginally

cor(final.df$flo_time,final.df$av_fruit_time)


###full phylogenetically, with hysteranthy to include synanthy
full.modAA<-phyloglm(pro2~pol+class2+shade_bin+fruit_bin+flo_type,final.df, pruned.by.anthy, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot = 0, full.matrix = TRUE)
summary(full.modAA)

###model with early subset only
earl<-  earl %>% remove_rownames %>% column_to_rownames(var="name")
earl.mod<-phyloglm(pro~pol+heigh_height+shade_bin+flo_type,earl, pruned.earl, method = "logistic_MPLE", btol = 10, log.alpha.bound = 4,
                    start.beta=NULL, start.alpha=NULL,
                    boot=100,full.matrix = TRUE)
summary(earl.mod)### yay, with onlu "early" fruiters

#########That was fun####################Nowdoit in BRMS############################################

library("brms")
library("MCMCglmm")
#https://cran.r-project.org/web/packages/brms/README.html for some guideance

#make all variable factoror brms
final.df$pro<-as.factor(final.df$pro)
final.df$pol<-as.factor(final.df$pol)
final.df$class2<-as.factor(final.df$class2)
final.df$shade_bin<-as.factor(final.df$shade_bin)
final.df$fruit_bin<-as.factor(final.df$fruit_bin)
final.df$flo_type<-as.factor(final.df$flo_type)

##construct covarience matrix:
final.df<-rownames_to_column(final.df, "name")
inv.phylo <- MCMCglmm::inverseA(pruned.by.anthy, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

###Best model###############################################################################
model <- brm(pro~ pol+class2+fruit_bin+shade_bin +flo_type+ (1|name), data = final.df, 
 family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
 prior = c(prior(normal(0, 5), "b"),
 prior(normal(0, 5), "Intercept"),
 prior(student_t(3, 0, 5), "sd"))) ###why does sigma not appear in my model?? #should list(name or pruned.by.anthy)
summary(model)

###bayesian and continuous
modelcont <- brm(pro~ pol+heigh_height+av_fruit_time+shade_bin +flo_type+ (1|name), data = final.df, 
             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
             prior = c(prior(normal(0, 5), "b"),
                       prior(normal(0, 5), "Intercept"),
                       prior(student_t(3, 0, 5), "sd"))) 
summary(full.modB)

####Phylogenetic signal:Work in progress
hyp<-"sd_name__Intercept^2/(sd_name__Intercept^2+((3.14159^2)/3)) =0" ###This might be the phylogenetic correlation: https://stats.stackexchange.com/questions/62770/calculating-icc-for-random-effects-logistic-regression
hypothesis(model,hyp, class=NULL)
#explainaition of iCC http://www.theanalysisfactor.com/the-intraclass-correlation-coefficient-in-mixed-models/

summary(full.modA)
plot(marginal_effects(model, probs = c(0.05, 0.95)))

########################Visualization####################
load("/Users/danielbuonaiuto/Desktop/shinystan-multiparam-gg (2).RData")
library(ggplot2)
p<-shinystan_multiparam_gg

my_labels<-c("sd Name (intercept","flower type","shade tolerance", "height class", "dispersal season","pollination syndrome","intercept")
p+scale_y_continuous(breaks= 1:7,labels = my_labels)+theme_classic()+geom_vline(aes(xintercept=0))+labs(x="effect size", y="predictor")
### View the stan code for the mdoe
stancode(modelcont)
##predictions, not totally usefule
####################check out the priors############################################
beta_draws <- as.matrix(modelcont, pars = "b")
dim(beta_draws)
library(bayesplot)
mcmc_intervals(beta_draws)
beta2_and_prior <- cbind(
  prior = rnorm(nrow(beta_draws), 0, 10), # draw from prior distribution
  posterior = beta_draws[, 2]
)
mcmc_areas(beta2_and_prior) 
###try it with other values
in_draws <- as.matrix(modelcont, pars = "Intercept")
dim(in_draws)
mcmc_intervals(in_draws)
beta3_and_prior <- cbind(
  prior = rnorm(nrow(in_draws), 0, 10), # draw from prior distribution
  posterior = beta_draws[, 2]
)
mcmc_areas(beta3_and_prior) 
###student T
sd_draws <- as.matrix(modelcont, pars = "sd")
dim(sd_draws)
mcmc_intervals(sd_draws)
beta4_and_prior <- cbind(
  prior = rnorm(nrow(sd_draws), 0, 10), # draw from prior distribution
  posterior = beta_draws[, 2]
)
mcmc_areas(beta4_and_prior)

#all together now

#mcmc_areas(beta2_and_prior)
#mcmc_areas(beta3_and_prior)
#mcmc_areas(beta4_and_prior)
### seems okay but ask Lizzie
########################posterior###check#######################################################
plot(model)
pp_check(model, type = "bars")


#####see in shiny################################
library(shinystan)                   
launch_shiny(model, rstudio = getOption("shinystan.rstudio"))
launch_shiny(modelcont, rstudio = getOption("shinystan.rstudio"))


stop("no need to run the full crossed model")
##########################full crossed model#############################################################
model2 <- brm(pro~ pol*class2*fruit_bin*shade_bin *flo_type+ (1|name), data = final.df, 
             family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=5000,
             prior = c(prior(normal(0, 5), "b"),
                       prior(normal(0, 5), "Intercept"),
                       prior(student_t(3, 0, 5), "sd")))
summary(model2) ### no interactions significant


table(anthy$Phen.sequence)
library(rstan)




