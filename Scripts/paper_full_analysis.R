####This is the master Rcode for Dan B's hysteranthy paper

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")

###libraryies
library("tidyverse")
library("ggplot2")
library(brms)
library(ggthemes)
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library(boot)
library("phylolm")
library(ggstance)

library("raster")
library("remote")
library(reshape2)
library(RColorBrewer)

load("RData/paper_full_analysis.RData")
 ###Data#####################################################################
#Read in Pep data
aln<-read.csv("datasheets_derived/alnus_delta_hyst.csv",header=TRUE)
frax<-read.csv("datasheets_derived/fraxinus_delta_hyst.csv",header=TRUE)
aes<-read.csv("datasheets_derived/aes_delta_hyst.csv",header=TRUE)

# Readin Harard forest data 2,3, and 4c
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
HF$phys.offset<-HF$bb.jd-HF$fbb.jd
HF$funct.offset<-HF$l75.jd-HF$fopn.jd

###Data 4a, 4b but outdated, read in below when analysis begins
#mich.tree<-read.tree("pruned_for_mich.tre")
#mich.data<-read.csv("mich_data_full_clean.csv")
#drought.dat<-read.csv("..//Data/USDA_traitfor_MTSV.csv",header=TRUE)

##Read in larger pep data
aln10<-read.csv("datasheets_derived/alnus10_delta_hyst.csv",header=TRUE)
frax10<-read.csv("datasheets_derived/fraxinus10_delta_hyst.csv",header=TRUE)
bet10<-read.csv("datasheets_derived/betpen10_delta_hyst.csv",head=TRUE)
aes10<-read.csv("datasheets_derived/aes10_delta_hyst.csv",header=TRUE)

aln10$taxa<-"Alnus glutinosa" ##assign species
frax10$taxa<-"Fraxinus excelsior"
aes10$taxa<-"Aesculus hippocastenum"

intra.d<-rbind(aln10,frax10,bet10,aes10)
moist.aug<-raster("grids_germany_multi_annual_soil_moist_1991-2010_08.asc") ##Gause Kruger 3 for August
moist.apr<-raster("grids_germany_multi_annual_soil_moist_1991-2010_04.asc")

###analysis 1: Is hysteranthy changing with time? featureing PEP725 data

###function for ordering pepsite
pepnumber <- function(dat, sitecolname){
  df <- data.frame(s_id=unique(as.numeric(unlist(dat[sitecolname]))),
                   peporder=c(1:length(unique(unlist(dat[sitecolname])))))
  datmerge <- merge(dat, df, by=sitecolname)
  return(datmerge)
  
}

###making the hinge
aln <- pepnumber(aln, "s_id")

aln$YEAR.hin <- aln$year
aln$YEAR.hin[which(aln$YEAR.hin<1980)] <- 1980
aln$YEAR.hin <- aln$YEAR.hin-1980

frax <- pepnumber(frax, "s_id")
frax$YEAR.hin <- frax$year
frax$YEAR.hin[which(frax$YEAR.hin<1980)] <- 1980
frax$YEAR.hin <- frax$YEAR.hin-1980

aes<- pepnumber(aes, "s_id")
aes$YEAR.hin <- aes$year
aes$YEAR.hin[which(aes$YEAR.hin<1980)] <- 1980
aes$YEAR.hin <- aes$YEAR.hin-1980

##### Alnus model
fit.aln.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aln) 
aln.dat<-as.data.frame(coef(fit.aln.brms))

aln.dat<-rownames_to_column(aln.dat, var = "peporder")
aln.dat<-merge(aln.dat,aln)
colnames(aln.dat)

names(aln.dat)[1]<-"peporder"
names(aln.dat)[4]<-"Intercept.2.5"
names(aln.dat)[5]<-"Intercept.97.5"
names(aln.dat)[2]<-"Intercept"
names(aln.dat)[6]<-"slope"
names(aln.dat)[7]<-"slope.2.5"
names(aln.dat)[8]<-"slope.97.5"
aln.dat<-dplyr::select(aln.dat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
aln.dat$hinge<-ifelse(aln.dat$YEAR.hin>0,1,0)

alphaALN<-mean(aln.dat$Intercept)
alphaALNlow<-mean(aln.dat$Intercept.2.5)
alphaALNhigh<-mean(aln.dat$Intercept.97.5)
betaALN<-mean(aln.dat$slope)
betaALNlow<-mean(aln.dat$slope.2.5)
betaALNhigh<-mean(aln.dat$slope.97.5)

####Fraxinus model
fit.frax.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=frax) 
frax.dat<-as.data.frame(coef(fit.frax.brms))
frax.dat<-rownames_to_column(frax.dat, var = "peporder")
frax.dat<-merge(frax.dat,frax)

names(frax.dat)[1]<-"peporder"
names(frax.dat)[4]<-"Intercept.2.5"
names(frax.dat)[5]<-"Intercept.97.5"
names(frax.dat)[2]<-"Intercept"
names(frax.dat)[6]<-"slope"
names(frax.dat)[7]<-"slope.2.5"
names(frax.dat)[8]<-"slope.97.5"
frax.dat<-dplyr::select(frax.dat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
frax.dat$hinge<-ifelse(frax.dat$YEAR.hin>0,1,0)


alphaFRAX<-mean(frax.dat$Intercept)
alphaFRAXlow<-mean(frax.dat$Intercept.2.5)
alphaFRAXhigh<-mean(frax.dat$Intercept.97.5)
betaFRAX<-mean(frax.dat$slope)


###Aesculus model
fit.aes.brms<-brm(offset~YEAR.hin+(YEAR.hin|peporder),data=aes) 
aes.dat<-as.data.frame(coef(fit.aes.brms))
aes.dat<-rownames_to_column(aes.dat, var = "peporder")
aes.dat<-merge(aes.dat,aes)

names(aes.dat)[1]<-"peporder"
names(aes.dat)[4]<-"Intercept.2.5"
names(aes.dat)[5]<-"Intercept.97.5"
names(aes.dat)[2]<-"Intercept"
names(aes.dat)[6]<-"slope"

names(aes.dat)[7]<-"slope.2.5"
names(aes.dat)[8]<-"slope.97.5"
aes.dat<-dplyr::select(aes.dat,peporder,Intercept,Intercept.2.5,Intercept.97.5,slope,slope.2.5,slope.97.5,s_id,year,offset,YEAR.hin,lat,lon)
aes.dat$hinge<-ifelse(aes.dat$YEAR.hin>0,1,0)


alphaAES<-mean(aes.dat$Intercept)
betaAES<-mean(aes.dat$slope)
alphaAESlow<-mean(aes.dat$Intercept.2.5)
alphaAEShigh<-mean(aes.dat$Intercept.97.5)



allbetas<-aln.dat$slope
allalphas<-aln.dat$Intercept
allbetas1<-frax.dat$slope
allalphas1<-frax.dat$Intercept
allbetas2<-aes.dat$slope
allalphas2<-aes.dat$Intercept

# graph

png("..//figure/FLS_climate_change.png",width = 5, height = 6, units = 'in', res = 300)
plot(c(1960,2015), c(-30,50), type = "n", xlab = "year", ylab = "Days between flowering and leafing", bty='l')
#rect(xleft=1960, ybottom=-30, xright=1980, ytop=50,col="ivory1" )
#rect(xleft=1980, ybottom=-30, xright=2015, ytop=50,col="ivory2")
segments(x0=1980,y0=allalphas,x1=2015,y1=allalphas+allbetas*35,col=rgb(0.1,0.1,0.1,alpha=0.1),lty="dotted", lwd=0.1 )
segments(x0=1960,y0=allalphas,x1=1980,y1=allalphas,col=rgb(.1,.1,.1,alpha=0.1),lty="dotted",lwd=0.1 )
segments(x0=1980,y0=allalphas1,x1=2015,y1=allalphas1+allbetas1*35,col=rgb(0,0,0,alpha=0.1),lty="dotted",lwd=0.1)
segments(x0=1960,y0=allalphas1,x1=1980,y1=allalphas1,col=rgb(0.1,0.1,0.1,alpha=0.1),lty="dotted",lwd=0.1)
segments(x0=1980,y0=allalphas2,x1=2015,y1=allalphas2+allbetas2*35,col=rgb(0.1,0.1,0.1,alpha=0.1),lty="dotted",lwd=0.1 )
segments(x0=1960,y0=allalphas2,x1=1980,y1=allalphas2,col=rgb(0.1,0.1,0.1,alpha=0.1),lty="dotted",lwd=0.1 )
segments(x0=1960,y0=alphaALN, x1=1980,y1=alphaALN,col="red",lwd=3)
segments(x0=1980,y0=alphaALN, x1=2015,y1=alphaALN+betaALN*35,col="red",lwd=3)
rect(xleft=1960,xright=2015,ybottom=alphaALNlow,ytop=alphaALNhigh,col=rgb(1,0,0,alpha=0.3),border=NA)
#segments(x0=1980, y0=alphaALN-betaALNlow,x1=2015,y1=alphaALN-betaALNlow+betaALN*35) was trying to add cis to slope but too small to see

#segments(x0=1960,y0=alphaALNlow, x1=2015,y1=alphaALNlow, lty=4, col="red",lwd=1)
#segments(x0=1960,y0=alphaALNhigh, x1=2015,y1=alphaALNhigh, lty=4, col="red",lwd=1)
segments(x0=1960,y0=alphaFRAX, x1=1980,y1=alphaFRAX,col="darkgoldenrod1",lwd=3)
segments(x0=1980,y0=alphaFRAX, x1=2015,y1=alphaFRAX+betaFRAX*35,col="darkgoldenrod1",lwd=3)
rect(xleft=1960,xright=2015,ybottom=alphaFRAXlow,ytop=alphaFRAXhigh,col=rgb(1,.8,0,alpha=0.3),border=NA)

#segments(x0=1960,y0=alphaFRAXlow, x1=2015,y1=alphaFRAXlow, lty=4,col="yellow",lwd=1)
#segments(x0=1960,y0=alphaFRAXhigh, x1=2015,y1=alphaFRAXhigh, lty=4,col="yellow",lwd=1)
segments(x0=1960,y0=alphaAES, x1=1980,y1=alphaAES,col="blue",lwd=3)
segments(x0=1980,y0=alphaAES, x1=2015,y1=alphaAES+betaAES*35,col="blue",lwd=3)
rect(xleft=1960,xright=2015,ybottom=alphaAESlow,ytop=alphaAEShigh,col=rgb(0,0,1,alpha=0.3),border=NA)

#segments(x0=1960,y0=alphaAESlow, x1=2015,y1=alphaAESlow, lty=4,col="blue",lwd=1)
#segments(x0=1960,y0=alphaAEShigh, x1=2015,y1=alphaAEShigh, lty=4,col="blue",lwd=1)
par(xpd=FALSE)
legend("top", legend=c("A. glutinosa", "F. excelsior","A. hippocastanum"),col=c("red", "darkgoldenrod1","blue"), lwd=2, cex=0.6, horiz=TRUE)
dev.off()

##################################Analysis II####################################################
unique(HF$species)
#rename all taxa
HF$name[HF$species=="ACPE"]<-"A. pensylvanicum"
HF$name[HF$species=="ACRU"]<-"*A. rubrum"
HF$name[HF$species=="ACSA"]<-"*A. saccharrum"
HF$name[HF$species=="AMSP"]<-"A. species"
HF$name[HF$species=="BEAL"]<-"*B. allegheniensis"
HF$name[HF$species=="BELE"]<-"*B. lenta"
HF$name[HF$species=="BEPA"]<-"*B. papyrifera"
HF$name[HF$species=="BEPO"]<-"*B. populifolia"
HF$name[HF$species=="FAGR"]<-"*F. grandifolia"
HF$name[HF$species=="FRAM"]<-"*F. americana"
HF$name[HF$species=="HAVI"]<-"H. virginiana"
HF$name[HF$species=="ILVE"]<-"I. verticillata"
HF$name[HF$species=="KAAN"]<-"K. angustifolia"
HF$name[HF$species=="KALA"]<-"K. latifolia"
HF$name[HF$species=="NEMU"]<-"I. mucronata"
HF$name[HF$species=="NYSY"]<-"N. sylvatica"
HF$name[HF$species=="POTR"]<-"*P. tremuloides"
HF$name[HF$species=="PRSU"]<-"P. serotina"
HF$name[HF$species=="QURU"]<-"*Q. rubra"
HF$name[HF$species=="QUVE"]<-"*Q. velutina*"


HF.averages<-dplyr::filter(HF, species %in% c("ACRU","BEAL" ,"QURU","ACPE","NYSY" ))

#HF.averages$name[HF.averages$species=="ACPE"]<-"A. pensylvanicum"
#HF.averages$name[HF.averages$species=="ACRU"]<-"A. rubrum"
#HF.averages$name[HF.averages$species=="BEAL"]<-"B. alleghaniensis"
#HF.averages$name[HF.averages$species=="NYSY"]<-"N. sylvatica"
#HF.averages$name[HF.averages$species=="QURU"]<-"Q. rubra"


colnames(HF.averages)
colnames(HF.averages)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion (75%)","flower budburst","flower open","physiological offset","functional offset","Funct. hysteranthy","Phys. hysteranthy","name")
HF.averages<-tidyr::gather(HF.averages, phase,DOY,4:7)

HF.averages$name <- factor(HF.averages$name, levels = c("A. rubrum","B. alleghaniensis" ,"Q. rubra","A. pensylvanicum","N. sylvatica"))
pd<-position_dodge(width=0.0)

###if you want a bigger data set to plot
colnames(HF)
colnames(HF)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion (75%)","flower budburst","flower open","physiological offset","functional offset","name")
HF2<-tidyr::gather(HF, phase,DOY,4:7)
HF2<-filter(HF2,!is.na(name))
png("..//figure/HFmeans_expanded.png",width = 8, height = 6, units = 'in', res=300)
ggplot(HF2,(aes(name,DOY)))+stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max,aes(shape=phase,color=phase))+scale_color_manual(values=c("darkgray","darkgray","black","black"))+scale_shape_manual(values=c(2,1,17,16))+theme_base()+ylab("Day of Year")+xlab(NULL)+theme(axis.text.x = element_text(angle = 300,hjust=0))
dev.off()
###reduce sp
#jpeg("..//figure/HFmeans.jpeg",width = 800, height = 550)
#ggplot(HF.averages,(aes(name,DOY)))+stat_summary(fun.data = "mean_cl_boot",aes(color=phase,shape=phase),position=pd,size=0.78)+scale_color_manual(values=c("firebrick1","deeppink","lawngreen","darkgreen"))+scale_shape_manual(values=c(16,8,0,23))+theme_base()+ylab("Day of Year")+xlab(NULL)#+theme(axis.text.x = element_text(angle = 300,hjust=0.5))
#dev.off()

########Analysis 3 Quercus rubra at HF.........
QURU<-filter(HF,species=="QURU")
QURU<-filter(QURU,tree.id!="QURU-02")
colnames(QURU)
QURU$FLS<-ifelse(QURU$"physiological offset"<0,"seranthous","hysteranthous")
colnames(QURU)<-c("year","tree.id","species","leaf_budburst","leaf_expansion(75%)","flower_budburst","flower_open","physiological_offset","functional_offset","name","FLS")
QURU<-drop_na(QURU)
QURU$tree.id[QURU$tree.id=="QURU-01"]<-"1"
QURU$tree.id[QURU$tree.id=="QURU-03"]<-"3"
QURU$tree.id[QURU$tree.id=="QURU-04"]<-"4"

jpeg("..//figure/HF_Q_ru_interannual.jpeg",width = 1500, height = 750,res=250)
pd2<-position_dodge2(0.8)

ggplot(QURU,aes(year,flower_budburst))+geom_point(aes(year,flower_budburst,color=tree.id,group = row.names(QURU)),shape=2, ,position=pd2)+
  geom_point(aes(year,leaf_budburst,color=tree.id,group = row.names(QURU)),shape=17,position=pd2)+geom_linerange(aes(x=year,ymin=flower_budburst,ymax=leaf_budburst, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd2)+theme_base()+labs(y = "Day of year",color= "Tree I.D.")+scale_color_manual(values=c("red","blue","darkgreen"))
dev.off()

#######Analysis 3 come back to this once Lizzie talk to Jonath.
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
  }
extract_CIs<-function(x){
  filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
  }  
###read in data
mich.data<-read.csv("datasheets_derived/MTSV_USFS/michdata_final.csv")
mich.tre<-read.tree("datasheets_derived/MTSV_USFS/michtre_final.tre")

silv.data<-read.csv("datasheets_derived/MTSV_USFS/silvdata_final.csv")
silv.tre<-read.tree("datasheets_derived/MTSV_USFS/silvtre_final.tre")
##### phylo signal for suppliment
mich.tre$node.label<-NULL
silv.tre$node.label<-NULL
d<-comparative.data(mich.tre,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
#physiological hysteranthy
Phylo.Pro.phys<-phylo.d(d, binvar=pro3) 
Phylo.Pro.phys
##functionalhysteranthy
PhyloPro.funct<-phylo.d(d,binvar=pro2)
PhyloPro.funct
plot(PhyloPro.funct)

d.silv<-comparative.data(silv.tre,silv.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
silv.Pro.phys<-phylo.d(d.silv, binvar=pro3) 
silv.Pro.phys
##functionalhysteranthy
silv.Pro.funct<-phylo.d(d.silv,binvar=pro2)
silv.Pro.funct
plot(PhyloPro.funct)

######### phyloglm requires species names to be in rownames
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<- silv.data %>% remove_rownames %>% column_to_rownames(var="name")

##without Interactions
z.funct.drought.noint<-phyloglm(pro2~pol_cent+flo_cent+precip_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)

z.phys.drought.noint<-phyloglm(pro3~pol_cent+flo_cent+precip_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

#With
#z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          #start.beta=NULL, start.alpha=NULL,
                          #boot=599,full.matrix = TRUE)


#z.phys.drought<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
 #                        start.beta=NULL, start.alpha=NULL,
  #                       boot=599,full.matrix = TRUE)

z.funct.drought.silvics.noint<-phyloglm(pro2~pol_cent+flo_cent+precip_cent,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                  start.beta=NULL, start.alpha=NULL,
                                  boot=599,full.matrix = TRUE)

z.phys.drought.silvics.noint<-phyloglm(pro3~pol_cent+flo_cent+precip_cent,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                 start.beta=NULL, start.alpha=NULL,
                                 boot=599,full.matrix = TRUE)
###side bar compare phylo to brms###########

#inv.phylo <- MCMCglmm::inverseA(mich.tre, nodes = "TIPS", scale = TRUE)
#A <- solve(inv.phylo$Ainv)
#rownames(A) <- rownames(inv.phylo$Ainv)

###

#mich.data$name<-rownames(mich.data)

#modelcont.funct.noint <- brm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent+(1|name), data = mich.data, 
#                             family = bernoulli(link = "logit"), cov_ranef = list(name= A),iter=3000) 

#binaryPGLMM(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre)

#z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
 #                         start.beta=NULL, start.alpha=NULL,
#                         boot=599,full.matrix = TRUE)
#summary( z.funct.drought)

##prep for plotting
mich.funct.noint.dat<-full_join(extract_coefs(z.funct.drought.noint),extract_CIs(z.funct.drought.noint),by="trait")
colnames(mich.funct.noint.dat)<-c("trait","estimate","low","high")
mich.funct.noint.dat$class<-"functional-MTSV"

mich.phys.noint.dat<-full_join(extract_coefs(z.phys.drought.noint),extract_CIs(z.phys.drought.noint),by="trait")
colnames(mich.phys.noint.dat)<-c("trait","estimate","low","high")
mich.phys.noint.dat$class<-"physiological-MTSV"

silv.funct.noint.dat<-full_join(extract_coefs(z.funct.drought.silvics.noint),extract_CIs(z.funct.drought.silvics.noint),by="trait")
colnames(silv.funct.noint.dat)<-c("trait","estimate","low","high")
silv.funct.noint.dat$class<-"functional-USFS"

silv.phys.noint.dat<-full_join(extract_coefs(z.phys.drought.silvics.noint),extract_CIs(z.phys.drought.silvics.noint),by="trait")
colnames(silv.phys.noint.dat)<-c("trait","estimate","low","high")
silv.phys.noint.dat$class<-"physiological-USFS"

michigan.noint<-rbind(mich.phys.noint.dat,mich.funct.noint.dat)
USFS.noint<-rbind(silv.phys.noint.dat,silv.funct.noint.dat)

comps<-rbind(michigan.noint,USFS.noint)
comps$category<-NA
comps$category[which(comps$class=="physiological-USFS")] <- "physiological"
comps$category[which(comps$class=="physiological-MTSV")] <- "physiological"
comps$category[which(comps$class=="functional-USFS")] <- "functional"
comps$category[which(comps$class=="functional-MTSV")] <- "functional"

comps$data<-NA
comps$data[which(comps$class=="physiological-USFS")] <- "USFS"
comps$data[which(comps$class=="physiological-MTSV")] <- "MTSV"
comps$data[which(comps$class=="functional-USFS")] <- "USFS"
comps$data[which(comps$class=="functional-MTSV")] <- "MTSV"

###change the variable names
comps$trait[which(comps$trait=="pol_cent")] <- "pollination syndrome"
comps$trait[which(comps$trait=="flo_cent")] <- "flowering time"
comps$trait[which(comps$trait=="precip_cent")] <- "minimum precipitation"


###looks better on seperate plots no it doesnt
#comps.MTSV<-filter(comps,data=="MTSV")
#comps.USFS<-filter(comps,data=="USFS")
pd=position_dodgev(height=0.4)


#jpeg("..//figure/MTSV_effectsize.jpeg",width = 800, height = 350)
#ggplot(comps.MTSV,aes(estimate,trait))+geom_point(size=4,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)
#dev.off()

#jpeg("..//figure/USFS_effectsize.jpeg",width = 800, height = 350)
#ggplot(comps.USFS,aes(estimate,trait))+geom_point(size=4,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_base(base_size=11)+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)
#dev.off()

##both



#
###Analysis 4cd
source("..//Scripts/continuous_mod_prep.R")
HF<-filter(HF,species!=("QUAL"))
spforcontmods<-df$species

HF.continuous.data<-filter(HF,species %in% c(spforcontmods))
HF.continuous.data<-left_join(HF.continuous.data,traits, by="species") ###This is the data for the continuous models

##zscore predictors for these models
HF.continuous.data$cent_pol<-(HF.continuous.data$pol-mean(HF.continuous.data$pol,na.rm=TRUE))/(2*sd(HF.continuous.data$pol,na.rm=TRUE))
HF.continuous.data$cent_minP<-(HF.continuous.data$min_precip-mean(HF.continuous.data$min_precip))/(2*sd(HF.continuous.data$min_precip))
HF.continuous.data$cent_floday<-(HF.continuous.data$fopn.jd-mean(HF.continuous.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.continuous.data$fopn.jd,na.rm=TRUE))

###phylo lm can only handle one value per species in the tree
#pol.sum<- HF.continuous.data %>% group_by(name) %>%summarise(mean_pol=mean(cent_pol,na.rm=TRUE))
#flo.sum<- HF.continuous.data %>% group_by(name)%>%summarise(mean_floday=mean(cent_floday,na.rm=TRUE))
#P.sum<- HF.continuous.data %>% group_by(name)%>%summarise(mean_minP=mean(cent_minP,na.rm=TRUE))
#funct.sum<-HF.continuous.data %>% group_by(name)%>%summarise(mean_funct_offset=mean(funct.offset,na.rm=TRUE))
#phys.sum<-HF.continuous.data %>% group_by(name)%>%summarise(mean_phys_offset=mean(phys.offset,na.rm=TRUE))
#mean.continuous.dat<-left_join(funct.sum,phys.sum)
#mean.continuous.dat<-left_join(mean.continuous.dat,P.sum)
#mean.continuous.dat<-left_join(mean.continuous.dat,flo.sum)
#mean.continuous.dat<-left_join(mean.continuous.dat,pol.sum)



#df.cont<-mean.continuous.dat[match(mytree.names, mean.continuous.dat$name),]
#setdiff(mean.continuous.dat$name,mytree.names)
#setdiff(mytree.names,df.cont$name)
#df.cont$name==mytree.names
#df.cont$funct.bin<-ifelse(df.cont$mean_funct_offset>=0,1,0)
#df.cont$phys.bin<-ifelse(df.cont$mean_phys_offset>=0,1,0)
#df.cont<-as.data.frame(df.cont)


###HF phylosig
d3<-comparative.data(HF.tree.pruned,df.cont,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
#physiological hysteranthy
#phylo.d(d3, binvar=phys.bin) 
#phylo.d(d3, binvar=funct.bin)
#df.cont<- df.cont%>% remove_rownames %>% column_to_rownames(var="name")


#write.csv(df.cont,"HarvardForest/HF.means.data.csv")
#write.tree(HF.tree.pruned,"HarvardForest/HFtree.tre")
### continuous models
#HF.tree.pruned$tip.label
#rownames(df.cont)
#mean.continuous.dat$name
#HF.cont.funct<-phylolm(mean_funct_offset~mean_pol+mean_floday+mean_minP,data=df.cont,phy=HF.tree.pruned,boot=599)
#summary(HF.cont.funct)

#HF.cont.phys<-phylolm(mean_phys_offset~mean_pol+mean_floday+mean_minP,data=df.cont,phy=HF.tree.pruned,boot=599)
#summary(HF.cont.phys)

#HF.fuct.continuous<-full_join(extract_coefs(HF.cont.funct),extract_CIs(HF.cont.funct),by="trait")
#colnames(HF.fuct.continuous)<-c("trait","estimate","low","high")
#HF.fuct.continuous$category<-"functional"

#HF.phys.continuous<-full_join(extract_coefs(HF.cont.phys),extract_CIs(HF.cont.phys),by="trait")
#colnames(HF.phys.continuous)<-c("trait","estimate","low","high")
#HF.phys.continuous$category<-"physiological"


#HF.cont.comps<-rbind(HF.phys.continuous,HF.fuct.continuous)
#HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_pol")] <- "pollination syndrome"
#HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_floday")] <- "flowering time"
#HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_minP")] <- "minimum precipitation"


#HF.cont.comps<-filter(HF.cont.comps,trait!="sigma2")

#jpeg("..//figure/HF_cont_effectsize.jpeg",width = 800, height = 350)
#ggplot(HF.cont.comps,aes(estimate,trait))+geom_point(size=4,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))
#dev.off()

#HF.cont.comps$data<-"HF-continuous"
#comps<-dplyr::select(comps,-class)

#threeway<-rbind(HF.cont.comps,comps)

#ggplot(threeway,aes(estimate,trait))+geom_point(size=4,aes(color=data,shape=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=data, linetype=category))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)

#HF.bin.funct<-phyloglm(funct.bin~mean_pol+mean_floday+mean_minP,df.cont, HF.tree.pruned, method = "logistic_MPLE", btol = 105, log.alpha.bound = 8,
  #                        start.beta=NULL, start.alpha=NULL,
   #                       boot=599,full.matrix = TRUE)
#summary(HF.bin.funct)

#HF.bin.phys<-phyloglm(phys.bin~mean_pol+mean_floday+mean_minP,df.cont, HF.tree.pruned, method = "logistic_MPLE", btol = 500, log.alpha.bound = 10,
 #                     start.beta=NULL, start.alpha=NULL,
    #                  boot=599,full.matrix = TRUE)
#summary(HF.bin.phys)

#HF.fuct.bin<-full_join(extract_coefs(HF.bin.funct),extract_CIs(HF.bin.funct),by="trait")
#colnames(HF.fuct.bin)<-c("trait","estimate","low","high")
#HF.fuct.bin$category<-"functional"

#HF.phys.bin<-full_join(extract_coefs(HF.bin.phys),extract_CIs(HF.bin.phys),by="trait")
#colnames(HF.phys.bin)<-c("trait","estimate","low","high")
#HF.phys.bin$category<-"physiological"


#HF.bin.comps<-rbind(HF.phys.bin,HF.fuct.bin)
#HF.bin.comps$data<-"Hf-binary"
#HF.bin.comps$trait[which(HF.bin.comps$trait=="mean_pol")] <- "pollination syndrome"
#HF.bin.comps$trait[which(HF.bin.comps$trait=="mean_floday")] <- "flowering time"
#HF.bin.comps$trait[which(HF.bin.comps$trait=="mean_minP")] <- "minimum precipitation"
  
#fourway<-rbind(threeway, HF.bin.comps)

#jpeg("..//figure/allcases.jpeg",width = 800, height = 350)
#ggplot(fourway,aes(estimate,trait))+geom_point(size=4,aes(color=data,shape=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=data,linetype=category))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)
#dev.off()



inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

###continuous models
modelcont.funct.noint <- brm(funct.offset~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                             family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelcont.phys.noint <- brm(phys.offset~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                            family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelbin.funct.noint <- brm(hyst.funct~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                            family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) ##1 divergent transition
#summary(modelbin.funct.noint)
modelbin.phys.noint <- brm(hyst.phys~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                           family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,.25,.75,0.975))),"trait")
}

funct.cont<-extract_coefs4HF(modelcont.funct.noint)
funct.cont$class<-"functional"
funct.bin<-extract_coefs4HF(modelbin.funct.noint)
funct.bin$class<-"functional"

phys.cont<-extract_coefs4HF(modelcont.phys.noint)
phys.bin<-extract_coefs4HF(modelbin.phys.noint)
phys.cont$class<-"physiological"
phys.bin$class<-"physiological"


cont.noint<-rbind(phys.cont,funct.cont)
bin.noint<-rbind(phys.bin,funct.bin)
cont.noint$data_type<-"continuous"
bin.noint$data_type<-"binary"

hfboth<-rbind(cont.noint,bin.noint)
hfboth$trait[which(hfboth$trait=="cent_pol")] <- "pollination syndrome"
hfboth$trait[which(hfboth$trait=="cent_floday")] <- "flowering time"
hfboth$trait[which(hfboth$trait=="cent_minP")] <- "minimum precipitation"
hfboth$trait[which(hfboth$trait=="Intercept")] <- "(Intercept)"
colnames(hfboth)[colnames(hfboth)=="class"] <- "category"
colnames(hfboth)[colnames(hfboth)=="Estimate"] <- "estimate"

comps2<-comps
colnames(comps2)<-c("trait","estimate","Q2.5","Q97.5","class","category","data")
comps2<-dplyr::select(comps2,-class)
hfboth<-dplyr::select(hfboth,-Q25)
hfboth<-dplyr::select(hfboth,-Q75)
hfboth<-dplyr::select(hfboth,-Est.Error)

colnames(hfboth)
colnames(comps2)
hfboth$data<-"HF"
comps2$data_type<-"binary"
comps2<-comps2[,c(1,2,3,4,5,7,6)]


alleffectos<-rbind(comps2,hfboth)
alleffectos<-filter(alleffectos,trait!="(Intercept)")
pd=position_dodgev(height=0.8)

jpeg("..//figure/allmods_effectsizes_combined.jpeg",width=1200,height=650,res=180)
ggplot(alleffectos,aes(estimate,trait))+geom_point(aes(shape=data_type,color=data,fill=category),position=pd,size=3,stroke=1.5)+scale_shape_manual(values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey"))+
 geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data,alpha=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("solid","solid"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-30,30)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))

dev.off()

 bothHFs<-ggplot(hfboth,aes(estimate,trait))+geom_point(aes(shape=category,color=data),size=3,position=pd2)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=category,color=data_type),position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+
  theme_base(base_size = 11)+scale_color_manual(values=c("blue", "springgreen4"))+
  xlim(-30,30)+scale_shape_manual(values=c(1,2))+scale_linetype_manual(values=c("solid","solid"))

binos<-ggplot(comps,aes(estimate,trait))+
  geom_point(size=4,aes(shape=category,color=data),position=pd2)+
  geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,linetype=category,color=data))+geom_vline(aes(xintercept=0))+
  theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)+scale_linetype_manual(values=c("solid","solid"))

jpeg("..//figure/cases_2pannel.jpeg",width = 700, height = 550)
grid.arrange(binos,bothHFs,nrow=2)
dev.off()

pd2=position_dodgev(height=0.4)
hf.cont.noint<-ggplot(cont.noint,aes(Estimate,trait))+geom_point(aes(color=class),position=pd2)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=class),linetype="solid",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("HF-Continuous")
hf.bin.noint<-ggplot(bin.noint,aes(Estimate,trait))+geom_point(aes(color=class),position=pd2)+geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=class),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=class),linetype="dotted",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("HF-Binary")

png("..//figure/HF_cont_and_bin.png",)
grid.arrange(hf.cont.noint,hf.bin.noint,nrow=1) ###50 and 80 CIS
dev.off()

########Analysis 5############################Intra anuual
coordinates(intra.d)<- ~lon + lat ###PEP site coordinates
proj4string(intra.d) <- CRS("+proj=longlat +datum=WGS84") ##establish the project

#   <- spTransform(d, CRS("+proj=tmerc +lat_0=47.2700 +lon_0=7.5000 +k=1 +x_0=3386564.9400+y_0=5237917.9109 +ellps=krass +units=m +no_defs")) ### this doesnt work
p <- spTransform(intra.d, CRS("+proj=tmerc +lat_0=50.625 +lon_0=9.84375, +k=1 +x_0=3559832.734474 +y_0=5610860.187573 +ellps=krass +units=m +no_defs")) ###convert our coordinate to Gausse Kruger


gaas<-coordinates(p) ### make these values coordinates in structure

gaas.cord<-as.data.frame(gaas)
pep.cord<-as.data.frame(intra.d)
colnames(gaas.cord)<-c("x","y")


intra.df<-cbind(pep.cord,gaas.cord) ### aff the next

###Extract soil moisture baluse
Soil<-extract(moist.aug, matrix(c(intra.df$x,intra.df$y), ncol = 2)) ### extract the soil moisture at the pep sites
intra.df$SM<-Soil

###center phenology
intra.df$flo.cent<-intra.df$flower-mean(intra.df$flower,na.rm=TRUE)
intra.df$leaf.cent<-intra.df$leaf-mean(intra.df$leaf,na.rm=TRUE)
intra.df$soil.cent<-intra.df$SM-mean(intra.df$SM,na.rm=TRUE)
intra.df$offset.cent<-intra.df$offset-mean(intra.df$offset,na.rm=TRUE)
table(intra.df$soil.cent)
####now run each species model's seperately
#df.intra.alnus<-filter(intra.df,taxa=="Alnus glutinosa")
#df.intra.frax<-filter(intra.df,taxa=="Fraxinus excelsior")
#df.intra.bet<-filter(intra.df,taxa=="Betula pendula")
#df.intra.aes<-filter(intra.df,taxa=="Aesculus hippocastenum")





### leafing or flowering driving variation
colnames(intra.df)
mean(df.intra.alnus$offset)
summary(lm(offset~flo.cent,data=df.intra.alnus))
summary(lm(offset~leaf.cent,data=df.intra.alnus))


summary(lm(offset~flo.cent,data=df.intra.frax))
summary(lm(offset~leaf.cent,data=df.intra.frax))
summary(lm(offset~flo.cent,data=df.intra.aes))
summary(lm(offset~leaf.cent,data=df.intra.aes))
summary(lm(offset~flo.cent,data=df.intra.bet))
summary(lm(offset~leaf.cent,data=df.intra.bet))

par(mfrow=c(1, 2) )
plot(c(-10,10), c(-80,80), type = "n", xlab = "Flowering day (deviation from mean)", ylab = "change in FLS offset", bty='l')
segments(x=0,y0=6.964847 ,x1=10,y1=6.964847-.67*10,col="darkgreen",lty="solid" )
segments(x=0,y0=6.964847 ,x1=-10,y1=6.964847+.67*10,col="darkgreen",lty="solid" )
points(x=df.intra.alnus$flo.cent,y=df.intra.alnus$offset,col="darkgreen",pch=".",size=0.4)

segments(x=0,y0=12.193582 ,x1=10,y1=12.193582-.53*10,col="red",lty="solid" )
segments(x=0,y0=12.193582 ,x1=-10,y1=12.193582+.53*10,col="red",lty="solid" )
segments(x=0,y0=-14.127707 ,x1=10,y1=-14.127707-.28*10,col="blue",lty="solid" )
segments(x=0,y0=-14.127707 ,x1=-10,y1=-14.127707+.28*10,col="blue",lty="solid" )

plot(c(-10,10), c(-40,40), type = "n", xlab = "Leafing day (deviation from mean)", ylab = "chang in FLS offset", bty='l')
segments(x=0,y0=32.810102 ,x1=10,y1=32.810102+.23*10,col="darkgreen",lty="solid" )
segments(x=0,y0=32.810102 ,x1=-10,y1=32.810102-.23*10,col="darkgreen",lty="solid" )
segments(x=0,y0=6.746710  ,x1=10,y1=6.746710 +.24*10,col="red",lty="solid" )
segments(x=0,y0=6.746710  ,x1=-10,y1=6.746710 -.24*10,col="red",lty="solid" )
segments(x=0,y0=-17.618607 ,x1=10,y1=-17.618607+.4*10,col="blue",lty="solid" )
segments(x=0,y0=-17.618607 ,x1=-10,y1=-17.618607-.4*10,col="blue",lty="solid" )

###3 sopil moisture basyiasn
smprior<-get_prior(offset~soil.cent+flo.cent,data=intra.df)
#m.bays<-brm(offset~soil.cent+(1|taxa),data=intra.df,prior=smprior)
smflo.bayes<-brm(offset~flo.cent+soil.cent,data=intra.df,prior=smprior)
summary(smflo.bayes)
pepbayes<-extract_coefs4HF(smflo.bayes)
colnames(alleffectos)
colnames(pepbayes)
pepbayes<-dplyr::select(pepbayes,-Est.Error,-Q75,-Q25)
pepbayes$category<-"other"
pepbayes$data_type<-"continuous"
pepbayes$data<-"PEP725"
colnames(pepbayes)<-c("trait"  ,   "estimate",  "Q2.5"    ,  "Q97.5"   ,  "category",  "data_type" ,"data"  )
pepbayes$trait[which(pepbayes$trait=="flo.cent")] <- "flowering time"
pepbayes$trait[which(pepbayes$trait=="soil.cent")] <- "minimum precipitation"
pepbayes<-filter(pepbayes,trait!="Intercept")
biggest<-rbind(alleffectos,pepbayes)
biggest$trait[which(biggest$trait=="minimum precipitation")] <- "water dynamics"

jpeg("..//figure/allmods_effectsizes_combined.jpeg",width=1400,height=900,res=180)
ggplot(biggest,aes(estimate,trait))+geom_point(aes(shape=data_type,color=data,fill=category),position=pd,size=3,stroke=1.5)+scale_shape_manual(values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey",other="beige"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data,alpha=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("solid","solid"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-30,30)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4","skyblue1"))+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",other="beige",physiological="gray"))))
dev.off()

library(lme4)

sm.lmer<-lmer(offset~soil.cent+(1|taxa),data=intra.df)
smflo.lmer<-lmer(offset~flo.cent+soil.cent+(1|taxa),data=intra.df)
summary(sm.lmer)
ranef(sm.lmer)
newdat.lmer = data.frame(taxa = intra.df$taxa,
                        soil.cent = intra.df$soil.cent,
                        flo.cent = intra.df$flo.cent)

newdat.lmer$pred.sm<-predict(sm.lmer, newdata = newdat.lmer,re.form=NA)
newdat.lmer$pred.smflo<-predict(smflo.lmer, newdata = newdat.lmer,re.form=NA)
newdat.lmer<-gather(newdat.lmer,model,offset,4:5)
newdat.lmer$model[which(newdat.lmer$model=="pred.sm")] <- "soil moisture only"
newdat.lmer$model[which(newdat.lmer$model=="pred.smflo")] <- "soil moisture & flowering"

jpeg("..//figure/SM_comp.jpeg")
ggplot(newdat.lmer, aes(x = soil.cent, y = offset,linetype=model))+geom_smooth(method=lm,color="black",size=.3)+theme_classic()+xlab("deviation from mean soil moisture (% plant usable water)")+ylab("FLS offset")
dev.off()



save.image("paper_full_analysis.RData")
