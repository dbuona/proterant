rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/sub_proj/writing")
library(xtable,quietly=TRUE)
library(caper,quietly=TRUE)
library(dplyr,quietly=TRUE)
#read in data
mich.data<-read.csv("..//MTSV_USFS/michdata_final.csv")
mich.tre<-read.tree("..//MTSV_USFS/michtre_final.tre")

silv.data<-read.csv("..///MTSV_USFS/silvdata_final.csv")
silv.tre<-read.tree("..//MTSV_USFS/silvtre_final.tre")

HF<-read.csv("..//HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
HFsubber<-read.csv("..//HarvardForest/HFdata4modeling.csv",header=TRUE)
HF.tree<-read.tree("..//HarvardForest/HFtree4modeling.tre")
###make fls measure
HF$phys.fls<-HF$bb.jd-HF$fbb.jd
HF$funct.fls<-HF$l75.jd-HF$fopn.jd
HF$inter.fls<-HF$bb.jd-HF$fopn.jd

###make catagorical FLS
HF$hyst.funct<-ifelse(HF$funct.fls>0,1,0)
HF$hyst.phys<-ifelse(HF$phys.fls>0,1,0)
HF$hyst.inter<-ifelse(HF$inter.fls>0,1,0)

### prune the tree
HF<-dplyr::filter(HF,species!=("QUAL")) ## quercus alba has no flowers
spforcontmods<-HFsubber$species ##subset of species good for this analysis

HF.data<-dplyr::filter(HF,species %in% c(spforcontmods))
HF.data<-dplyr::left_join(HF.data,HFsubber, by="species") ###This is the data for the continuous models
meanhf<-HF.data %>% group_by(name) %>% summarise(meanfunctFLS=mean(funct.fls,na.rm=TRUE))
meanhf$FLSmeanfunctbin<-ifelse(meanhf$meanfunctFLS>0,1,0)
meanhf<-as.data.frame(meanhf)

HF.tree$node.label<-NULL
mich.tre$node.label<-NULL
silv.tre$node.label<-NULL


d.harv<-comparative.data(HF.tree,meanhf,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
Phylo.D.hf<-phylo.d(d.harv, binvar=FLSmeanfunctbin)

d<-comparative.data(mich.tre,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
Phylo.Pro.phys<-phylo.d(d, binvar=pro3) 

##functionalhysteranthy
PhyloPro.funct<-phylo.d(d,binvar=pro2)



d.silv<-comparative.data(silv.tre,silv.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
silv.Pro.phys<-phylo.d(d.silv, binvar=pro3) 
silv.Pro.funct<-phylo.d(d.silv,binvar=pro2)


##zscore predictors for these models
HF.data$pol_cent<-(HF.data$pol-mean(HF.data$pol,na.rm=TRUE))/(2*sd(HF.data$pol,na.rm=TRUE))
HF.data$precip_cent<-(HF.data$min_precip-mean(HF.data$min_precip))/(2*sd(HF.data$min_precip))
HF.data$flo_cent<-(HF.data$fopn.jd-mean(HF.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.data$fopn.jd,na.rm=TRUE))

HF.data$flo_cent.neg<--(HF.data$flo_cent)
HF.data$precip_cent.neg<--(HF.data$precip_cent)
###group by phylogeny
inv.phylo <- MCMCglmm::inverseA(HF.tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)


###continuous models normal
modelcont.funct <- brm(funct.fls~ pol+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol+pol:flo_cent+(1|name), data = HF.data, 
                       family = gaussian(), cov_ranef = list(name= A),iter=4000, warmup=3000) 

hyp <- "sd_name__Intercept^2 / (sd_name__Intercept^2 + sigma^2) = 0"

lambda.FLS <- hypothesis(modelcont.funct, hyp, class = NULL)


d2 <- density(lambda.FLS$samples[,1])


par(mfrow=c(3,2))
plot(Phylo.Pro.phys,main="MTSV-physiological")
plot(PhyloPro.funct,main="MTSV-functional")
plot(silv.Pro.phys, main="USFS-physiological")
plot(silv.Pro.funct,main="USFS-functional")
plot(Phylo.D.hf,main="HF-categorical")
plot(d2,main="HF-quantitative",xlab="Lambda",
     xlim=c(0,1),col="darkgray",fill="lightgray")
polygon(d2,col=adjustcolor("gray",0.4), border="gray")
abline(v=mean(lambda.FLS$samples[,1]),lty=1,col="black")

