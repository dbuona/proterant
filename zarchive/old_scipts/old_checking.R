### This is all old model checking code from concept.paper.analysis.R
###phyosig summary for 


(PhyloPro2) #0.06
PhyloPro3 #0.29
PhyloSilv2 #0.64
PhyloSilv3 #0.10


tab<-matrix(c(0.06,0.64,0.29,0.10))
rownames(tab)<-c("functional-MTSV","functional-Silvics","physiological-MTSV","Physiological-Silvics")
tab<-as.data.frame(tab)
tab<-rownames_to_column(tab, "Model")
###make this into a table at some point
###Which species are super early################################
superearl<-filter(mich.data, flo_time=="3.5")
superearl2<-filter(mich.data, flo_time=="4")
superearl3<-filter(mich.data, flo_time=="4.5")










stop("not an error, just ending the sourcing")


###########################################################################################

###average predictive comparsons(I dont think these are working yet.) Can I use divide by 4 rule?
coef(cent.funct.seed)
beta<-coef(cent.funct.seed)
hi<-1
lo<-0
###for [pollination syndrome]
beta[2]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$height_cent+beta[4]*mich.data$dev_time_cent+beta[5]*mich.data$fruit_cent+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$height_cent+beta[4]*mich.data$flo_cent+beta[5]*mich.data$dev_time_cent+beta[6]*mich.data$shade_bin)

print(mean(delta))  #0.41


###for flowering
april<-(4-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
may<-(5-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
earl<-2*sd(mich.data$flo_time)*april+mean(mich.data$flo_time)
mid<-2*sd(mich.data$flo_time)*may+mean(mich.data$flo_time)   

delta2<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height_cent+beta[4]*april+beta[5]*mich.data$dev_time_cent+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height_cent+beta[4]*may+beta[5]*mich.data$dev_time_cent+beta[6]*mich.data$shade_bin)

mean(delta2) ##.24 which is right! yay

2*sd(mich.data$flo_time)*(-0.697978)+mean(mich.data$flo_time)
2*sd(mich.data$flo_time)*(-0.1677311)+mean(mich.data$flo_time)                                            

###how do I unscale this to make it meaningful this?
#=#z score = X-m/sd
#rescale SD(z) + M
#I did mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
#so I should do 2sd(z)+M


###APC unscaled#########################################

summary(Mich.funct)
coef(Mich.funct)
beta<-coef(Mich.funct)
hi<-1
lo<-0
###for [pollination syndrome]
beta[2]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$heigh_height+beta[4]*mich.data$flo_time+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$heigh_height+beta[4]*mich.data$flo_time+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)

mean(delta)  #0.40

###for flowering 
earl<-4
mid<-5
delta2<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$heigh_height+beta[4]*earl+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$heigh_height+beta[4]*mid+beta[5]*mich.data$dev.time+beta[6]*mich.data$shade_bin)

mean(delta2)

#Whixh species are insect pollinated but not hysteranthous
buggy<-filter(mich.data, pol==0 & pro3==1)
unique(buggy$Family)
###tropical Rhamnaceae, Lauraceae, Anacardiaceae, Rutaceae
phy<-as.data.frame(coef(cent.funct.seed.winter))
phy<- phy %>%  rownames_to_column(var="effect")

nophy<-as.data.frame(coef(cent.funct.seed.nophylo))
nophy<- nophy %>%  rownames_to_column(var="effect")
phycom<-left_join(phy,nophy)

#-------------------------------------------------------#
#okay no to see if I can replicate the results of the scoop
###different hysteranthy definition
####different predictors
###didn't z-score (they log transformed continuous variables)
###reduced species list
#no interaactions
### they conflate early flowering hypotheses with hysteranthous hypotheses

### also I think they maybe ran a seperate model for each hypothesis

#1 ###take out flower time and compare functional to physiological
#predictor seed dev should be big and pol samall
funct.check<-phyloglm(pro2~pol+heigh_height+dev.time+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=10,full.matrix = TRUE)
summary(funct.check)
phys.check<-phyloglm(pro3~pol+heigh_height+dev.time+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                     start.beta=NULL, start.alpha=NULL,
                     boot=10,full.matrix = TRUE)
summary(phys.check)
##its true that dev time is increasing in importance and syndrome decreasing with different definitions
#### try it z scored
z.check.funct<-phyloglm(pro2~pol+height_cent++dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=10,full.matrix = TRUE)
summary(z.check.funct)
#pol 2.88, dev -0.5
z.check.phys<-phyloglm(pro3~pol+height_cent++dev_time_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=10,full.matrix = TRUE)
summary(z.check.phys)
#pol 2.01 dev -1.7

####lets look at dev time flower time and polinaton
z.check.funct2<-phyloglm(pro2~pol+flo_cent+dev_time_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=10,full.matrix = TRUE)

z.check.phys2<-phyloglm(pro3~pol+flo_cent+dev_time_cent,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                        start.beta=NULL, start.alpha=NULL,
                        boot=10,full.matrix = TRUE)
summary(z.check.funct2)
##pol 2.69, flo, -2.98, dev time -0.9
summary(z.check.phys2) 
##pol 1.5, flo -4.11, dev -1.43

funct.check2<-phyloglm(pro2~pol+flo_time+dev.time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                       start.beta=NULL, start.alpha=NULL,
                       boot=10,full.matrix = TRUE)
summary(funct.check2)
#pol 2.7, flo -1.86 dev -0.18

phys.check2<-phyloglm(pro3~pol+flo_time+dev.time,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=10,full.matrix = TRUE)
summary(phys.check2)
#pol 0.99, flo -2.46 dev -.22

####summary accross the board, polination syndrome decrease with phys definition

######do fruit and flo time covary? If yes removing one from the model should greatly change the effect
#test1<-phyloglm(pro2~pol+height_cent+flo_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#                         start.beta=NULL, start.alpha=NULL,
#                        boot=599,full.matrix = TRUE)
#summary(test1)

#test2<-phyloglm(pro2~pol+height_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
#               start.beta=NULL, start.alpha=NULL,
#              boot=599,full.matrix = TRUE)

#summary(test2)
####The effect sizes remain relatively stable there fore we assume there is not a collinearity issue
####marginal effects of flowering time on pollination

vels<-seq(0,1,1)
slopes <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[7]*vels
cbind(vels,slopes)
plot(vels, slopes, type = "l", lty = 1, ylim = c(-5, 0), xlab = "syndrome", ylab = "Marginal Effect of flower time")

###height
z.funct.seed.winter$coefficients
range(mich.data$height_cent)
height<-seq(-0.5603846, 1.0484774,0.1)
slopes2 <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[9]*height
cbind(height,slopes2)
plot(height, slopes2, type = "l", lty = 1, ylim = c(-7, 0), xlab = "height centered", ylab = "Marginal Effect of flower time")


###shade
z.funct.seed.winter$coefficients

shade<-seq(0,1,1)
slopes3 <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[10]*shade
cbind(shade,slopes3)
plot(shade, slopes3, type = "l", lty = 1, ylim = c(-4, 0), xlab = "shade centered", ylab = "Marginal Effect of flower time")

###seed development
z.funct.seed.winter$coefficients
range(mich.data$dev_time_cent)
seed<-seq(-0.7,2.1,0.2)
slopes4 <- z.funct.seed.winter$coefficients[4] + z.funct.seed.winter$coefficients[8]*seed
cbind(seed,slopes4)
plot(seed, slopes4, type = "l", lty = 1, ylim = c(-7, 2), xlab = "seed centered", ylab = "Marginal Effect of flower time")


###More plots for understanding interactions
interaction.plot(mich.data$pol,mich.data$flo_time,mich.data$pro2)
interaction.plot(mich.data$shade_bin,mich.data$flo_time,mich.data$pro2)

interaction.plot(mich.data$class,mich.data$flo_time,mich.data$pro2)

qplot(x = pol, y = pro2, data = mich.data, color = as.factor(flo_time)) +
  geom_smooth(method = "lm",se=FALSE), se=FALSE)+ggtitle("flower time X syndrome") 

qplot(x = shade_bin, y = pro2, data = mich.data, color = as.factor(flo_time)) +
  geom_smooth(method = "lm",se=FALSE)+ggtitle("flower time X shade tolerance") 

qplot(x = dev.time, y = pro2, data = mich.data, color = as.factor(flo_time)) +
  geom_smooth(method = "lm",se=FALSE)+ggtitle("flower time X development time") 

qplot(x = heigh_height ,y = pro2, data = mich.data, color = as.factor(flo_time)) +geom_smooth(method = "lm",se=FALSE)+ggtitle("flower time X height") 

###residuals and stuff

resid1<-as.data.frame(z.funct.seed.winter$fitted.values)
resid1<-rownames_to_column(resid1, "name")

resid2<-as.data.frame(z.funct.seed.winter$residuals)

resid2<-rownames_to_column(resid2, "name")
resid2$residuals<-resid$z.funct.seed.winter$residual
resid<-left_join(resid1,resid2)
qplot(z.funct.seed.winter$fitted.values,z.funct.seed.winter$residuals,data=resid)+geom_point()

goober<-filter(resid2,z.funct.seed.winter$residuals>=0.7)
goober2<-filter(resid2,z.funct.seed.winter$residuals<=-0.7)

residual.list<-rbind(goober,goober2)
write.csv(residual.list,"extreme.resid.csv")