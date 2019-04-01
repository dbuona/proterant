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

load("RData/paper_full_analysis.RData")
 ###Data#####################################################################
#data1
aln<-read.csv("datasheets_derived/alnus_delta_hyst.csv",header=TRUE)
frax<-read.csv("datasheets_derived/fraxinus_delta_hyst.csv",header=TRUE)
aes<-read.csv("datasheets_derived/aes_delta_hyst.csv",header=TRUE)

#data 2,3, and 4c
HF<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
HF$phys.offset<-HF$bb.jd-HF$fbb.jd
HF$funct.offset<-HF$l75.jd-HF$fopn.jd

###Data 4a, 4b
mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full_clean.csv")
drought.dat<-read.csv("..//Data/USDA_traitfor_MTSV.csv",header=TRUE)

##data5
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
#betaALNlow<-mean(aln.dat$slope.2.5)
#betaALNhigh<-mean(aln.dat$slope.97.5)

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

jpeg("..//figure/FLS_climate_change.jpeg")
plot(c(1960,2015), c(-30,45), type = "n", xlab = "year", ylab = "FLS offset", bty='l')
#rect(xleft=1960, ybottom=-30, xright=1980, ytop=50,col="ivory1" )
#rect(xleft=1980, ybottom=-30, xright=2015, ytop=50,col="ivory2")
segments(x0=1980,y0=allalphas,x1=2015,y1=allalphas+allbetas*35,col="azure3",lty="dotted" )
segments(x0=1960,y0=allalphas,x1=1980,y1=allalphas,col="azure3",lty="dotted")
segments(x0=1980,y0=allalphas1,x1=2015,y1=allalphas1+allbetas1*35,col="azure3")
segments(x0=1960,y0=allalphas1,x1=1980,y1=allalphas1,col="azure3")
segments(x0=1980,y0=allalphas2,x1=2015,y1=allalphas2+allbetas2*35,col="azure3",lty="dashed")
segments(x0=1960,y0=allalphas2,x1=1980,y1=allalphas2,col="azure3",lty="dashed")
segments(x0=1960,y0=alphaALN, x1=1980,y1=alphaALN,col="darkgreen",lwd=3)
segments(x0=1980,y0=alphaALN, x1=2015,y1=alphaALN+betaALN*35,col="darkgreen",lwd=3)
segments(x0=1960,y0=alphaALNlow, x1=2015,y1=alphaALNlow, lty=2, col="darkgreen",lwd=2)
segments(x0=1960,y0=alphaALNhigh, x1=2015,y1=alphaALNhigh, lty=2, col="darkgreen",lwd=2)
segments(x0=1960,y0=alphaFRAX, x1=1980,y1=alphaFRAX,col="red",lwd=3)
segments(x0=1980,y0=alphaFRAX, x1=2015,y1=alphaFRAX+betaFRAX*35,col="red",lwd=3)
segments(x0=1960,y0=alphaFRAXlow, x1=2015,y1=alphaFRAXlow, lty=2,col="red",lwd=2)
segments(x0=1960,y0=alphaFRAXhigh, x1=2015,y1=alphaFRAXhigh, lty=2,col="red",lwd=2)
segments(x0=1960,y0=alphaAES, x1=1980,y1=alphaAES,col="blue",lwd=3)
segments(x0=1980,y0=alphaAES, x1=2015,y1=alphaAES+betaAES*35,col="blue",lwd=3)
segments(x0=1960,y0=alphaAESlow, x1=2015,y1=alphaAESlow, lty=2,col="blue",lwd=2)
segments(x0=1960,y0=alphaAEShigh, x1=2015,y1=alphaAEShigh, lty=2,col="blue",lwd=2)
legend(2006,3, legend=c("A. glutinosa", "F. excelsior","A. hippocastanum"),col=c("darkgreen", "red","blue"), lwd=2, cex=0.6)
dev.off()

##################################Analysis II####################################################
HF.averages<-dplyr::filter(HF, species %in% c("ACRU","BEAL" ,"QURU","ACPE","NYSY" ))

HF.averages$name[HF.averages$species=="ACPE"]<-"A. pensylvanicum"
HF.averages$name[HF.averages$species=="ACRU"]<-"A. rubrum"
HF.averages$name[HF.averages$species=="BEAL"]<-"B. alleghaniensis"
HF.averages$name[HF.averages$species=="NYSY"]<-"N. sylvatica"
HF.averages$name[HF.averages$species=="QURU"]<-"Q. rubra"


colnames(HF.averages)
colnames(HF.averages)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion (75%)","flower budburst","flower open","physiological offset","functional offset","name")
HF.averages<-tidyr::gather(HF.averages, phase,DOY,4:7)

HF.averages$name <- factor(HF.averages$name, levels = c("A. rubrum","B. alleghaniensis" ,"Q. rubra","A. pensylvanicum","N. sylvatica"))
pd<-position_dodge(width=0.0)

jpeg("..//figure/HFmeans.jpeg")
ggplot(HF.averages,(aes(name,DOY)))+stat_summary(fun.data = "mean_cl_boot",aes(color=phase,shape=phase),position=pd,size=0.78)+scale_color_manual(values=c("deeppink","firebrick1","lawngreen","darkgreen"))+scale_shape_manual(values=c(0,8,16,23))+theme_base()+ylab("Day of Year")+xlab(NULL)#+theme(axis.text.x = element_text(angle = 300,hjust=0.5))
dev.off()

########Analysis 3 Quercus rubra at HF.........
QURU<-filter(HF,species=="QURU")
QURU<-filter(QURU,tree.id!="QURU-02")
QURU$FLS<-ifelse(QURU$phys.offset<0,"seranthous","hysteranthous")
QURU<-drop_na(QURU)
pd2<-position_dodge(0.6)

jpeg("..//figure/HF_Q.ru_interannual.jpeg")
ggplot(QURU,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd,color=tree.id,group = row.names(QURU)),shape=8, size=4,position=pd2)+geom_point(aes(year,bb.jd,color=tree.id,group = row.names(QURU)),shape=18,size=3,position=pd2)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd2)+theme_base()+labs(y = "Day of year",color= "Tree I.D.")
dev.off()

#######Analysis 3 come back to this once Lizzie talk to Jonath.
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait") ##This function extracts coefficients from phylolm model
  }
extract_CIs<-function(x){
  filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha") ##This function extracts CIs from phylo lm models.
  }  
###Quick clean
mich.data$pol<-ifelse(mich.data$Species=="quadrangulata",1,mich.data$pol)
mich.data$pol<-ifelse(mich.data$Genus=="Populus"& mich.data$Species=="nigra",1,mich.data$pol)
mich.tree$node.label<-NULL ###remove node labels so three works


###Analysis 4cd
source("..//Scripts/continuous_mod_prep.R")

spforcontmods<-traits$species

HF.continuous.data<-filter(HF,species %in% c(spforcontmods))
HF.continuous.data<-left_join(HF.continuous.data,traits, by="species") ###This is the data for the continuous models

##zscore predictors for these models
HF.continuous.data$cent_pol<-(HF.continuous.data$pol-mean(HF.continuous.data$pol,na.rm=TRUE))/(2*sd(HF.continuous.data$pol,na.rm=TRUE))
HF.continuous.data$cent_minP<-(HF.continuous.data$min_precip-mean(HF.continuous.data$min_precip))/(2*sd(HF.continuous.data$min_precip))
HF.continuous.data$cent_floday<-(HF.continuous.data$fopn.jd-mean(HF.continuous.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.continuous.data$fopn.jd,na.rm=TRUE))

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

modelbin.phys.noint <- brm(hyst.phys~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                           family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.10,.25,.75,0.90))),"trait")
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


pd2=position_dodgev(height=0.4)
hf.cont.noint<-ggplot(cont.noint,aes(Estimate,trait))+geom_point(aes(color=class),position=pd2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=class),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=class),linetype="dotted",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("HF-Continuous")
hf.bin.noint<-ggplot(bin.noint,aes(Estimate,trait))+geom_point(aes(color=class),position=pd2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=class),linetype="solid",position=pd2,width=0)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=class),linetype="dotted",position=pd2,width=0)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("orchid4", "springgreen4"))+ggtitle("HF-Binary")


jpeg("..//figure/HF_cont_and_bin.jpeg")
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

####now run each species model's seperately
df.intra.alnus<-filter(intra.df,taxa=="Alnus glutinosa")
df.intra.frax<-filter(intra.df,taxa=="Fraxinus excelsior")
df.intra.bet<-filter(intra.df,taxa=="Betula pendula")
df.intra.aes<-filter(intra.df,taxa=="Aesculus hippocastenum")

##model with just soil moisture
mod.intra.aln.sm<-brm(offset~soil.cent,data=df.intra.alnus)
summary(mod.intra.aln.sm)
mod.intra.frax.sm<-brm(offset~soil.cent,data=df.intra.frax)
summary(mod.intra.frax.sm)
mod.intra.bet.sm<-brm(offset~soil.cent,data=df.intra.bet)
summary(mod.intra.bet.sm)
mod.intra.aes.sm<-brm(offset~soil.cent,data=df.intra.aes)
summary(mod.intra.bet.sm)

priorz<-get_prior(offset~flo.cent*soil.cent,data=df.intra.alnus)
mod.intra.flo.sm.aln<-brm(offset~flo.cent*soil.cent,data=df.intra.alnus,prior=priorz)
summary(mod.intra.flo.sm.aln)
mod.intra.flo.sm.frax<-brm(offset~flo.cent*soil.cent,data=df.intra.frax,prior=priorz)
summary(mod.intra.flo.sm.frax)
mod.intra.flo.sm.sm<-brm(offset~flo.cent*soil.cent,data=df.intra.bet,prior=priorz)
summary(mod.intra.flo.sm.sm)

alnus.itra<-extract_coefs4HF(mod.intra.flo.sm.aln)
frax.itra<-extract_coefs4HF(mod.intra.flo.sm.frax)
bet.itra<-extract_coefs4HF(mod.intra.flo.sm.sm)
save.image("paper_full_analysis.RData")
