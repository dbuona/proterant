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
legend(2000,3, legend=c("A. glutinosa", "F. excelsior","A. hippocastanum"),col=c("darkgreen", "red","blue"), lwd=2, cex=0.6)
dev.off()

##################################Analysis II####################################################
HF.averages<-dplyr::filter(HF, species %in% c("ACRU","BEAL" ,"QURU","ACPE","NYSY" ))

HF.averages$name[HF.averages$species=="ACPE"]<-"A. pensylvanicum"
HF.averages$name[HF.averages$species=="ACRU"]<-"A. rubrum"
HF.averages$name[HF.averages$species=="BEAL"]<-"B. alleghaniensis"
HF.averages$name[HF.averages$species=="NYSY"]<-"N. sylvatica"
HF.averages$name[HF.averages$species=="QURU"]<-"Q. rubra"


colnames(HF.averages)
colnames(HF.averages)<-c("year" , "tree.id", "species","leaf budburst","leaf expansion (75%)","flower budburst","flower open","physiological offset","functional offset","Funct. hysteranthy","Phys. hysteranthy","name")
HF.averages<-tidyr::gather(HF.averages, phase,DOY,4:7)

HF.averages$name <- factor(HF.averages$name, levels = c("A. rubrum","B. alleghaniensis" ,"Q. rubra","A. pensylvanicum","N. sylvatica"))
pd<-position_dodge(width=0.0)

jpeg("..//figure/HFmeans.jpeg",width = 800, height = 550)
ggplot(HF.averages,(aes(name,DOY)))+stat_summary(fun.data = "mean_cl_boot",aes(color=phase,shape=phase),position=pd,size=0.78)+scale_color_manual(values=c("firebrick1","deeppink","lawngreen","darkgreen"))+scale_shape_manual(values=c(16,8,0,23))+theme_base()+ylab("Day of Year")+xlab(NULL)#+theme(axis.text.x = element_text(angle = 300,hjust=0.5))
dev.off()

########Analysis 3 Quercus rubra at HF.........
QURU<-filter(HF,species=="QURU")
QURU<-filter(QURU,tree.id!="QURU-02")
QURU$FLS<-ifelse(QURU$phys.offset<0,"seranthous","hysteranthous")
QURU<-drop_na(QURU)
QURU$tree.id[QURU$tree.id=="QURU-01"]<-"1"
QURU$tree.id[QURU$tree.id=="QURU-03"]<-"3"
QURU$tree.id[QURU$tree.id=="QURU-04"]<-"4"

jpeg("..//figure/HF_Q_ru_interannual.jpeg",width = 800, height = 550)
pd2<-position_dodge2(0.8)
ggplot(QURU,aes(year,fbb.jd))+geom_point(aes(year,fbb.jd,color=tree.id,group = row.names(QURU)),shape=8, size=4,position=pd2)+geom_point(aes(year,bb.jd,color=tree.id,group = row.names(QURU)),shape=18,size=3,position=pd2)+geom_linerange(aes(x=year,ymin=fbb.jd,ymax=bb.jd, linetype=FLS,color=tree.id,group = row.names(QURU)),position=pd2)+theme_base()+labs(y = "Day of year",color= "Tree I.D.")+scale_color_manual(values=c("red","blue","darkgreen"))
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

##with Interactions
z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)


z.phys.drought<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

z.funct.drought.silvics<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                  start.beta=NULL, start.alpha=NULL,
                                  boot=599,full.matrix = TRUE)

z.phys.drought.silvics<-phyloglm(pro3~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                 start.beta=NULL, start.alpha=NULL,
                                 boot=599,full.matrix = TRUE)
###side bar compare phylo to brms###########

inv.phylo <- MCMCglmm::inverseA(mich.tre, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

###

mich.data$name<-rownames(mich.data)

modelcont.funct.noint <- brm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent+(1|name), data = mich.data, 
                             family = bernoulli(link = "logit"), cov_ranef = list(name= A),iter=3000) 

binaryPGLMM(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre)

z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent+precip_cent+precip_cent:flo_cent+precip_cent:pol_cent+pol_cent:flo_cent,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)
summary( z.funct.drought)

##prep for plotting
mich.funct.wint.dat<-full_join(extract_coefs(z.funct.drought),extract_CIs(z.funct.drought),by="trait")
colnames(mich.funct.wint.dat)<-c("trait","estimate","low","high")
mich.funct.wint.dat$class<-"functional-MTSV"

mich.phys.wint.dat<-full_join(extract_coefs(z.phys.drought),extract_CIs(z.phys.drought),by="trait")
colnames(mich.phys.wint.dat)<-c("trait","estimate","low","high")
mich.phys.wint.dat$class<-"physiological-MTSV"

silv.funct.wint.dat<-full_join(extract_coefs(z.funct.drought.silvics),extract_CIs(z.funct.drought.silvics),by="trait")
colnames(silv.funct.wint.dat)<-c("trait","estimate","low","high")
silv.funct.wint.dat$class<-"functional-USFS"

silv.phys.wint.dat<-full_join(extract_coefs(z.phys.drought.silvics),extract_CIs(z.phys.drought.silvics),by="trait")
colnames(silv.phys.wint.dat)<-c("trait","estimate","low","high")
silv.phys.wint.dat$class<-"physiological-USFS"

michigan.wint<-rbind(mich.phys.wint.dat,mich.funct.wint.dat)
USFS.wint<-rbind(silv.phys.wint.dat,silv.funct.wint.dat)

comps<-rbind(michigan.wint,USFS.wint)
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
comps$trait[which(comps$trait=="pol_cent")] <- "main effect: pollination syndrome"
comps$trait[which(comps$trait=="flo_cent")] <- "main effect: flowering time"
comps$trait[which(comps$trait=="precip_cent")] <- "main effect: minimum precipitation"
comps$trait[which(comps$trait=="flo_cent:precip_cent")] <- "interaction: flowering x precip."
comps$trait[which(comps$trait=="pol_cent:precip_cent")] <- "interaction: pollination x precip."
comps$trait[which(comps$trait=="pol_cent:flo_cent")] <- "interaction: pollination x flowering"

###looks better on seperate plots
comps.MTSV<-filter(comps,data=="MTSV")
comps.USFS<-filter(comps,data=="USFS")
pd=position_dodgev(height=0.4)

jpeg("..//figure/MTSV_effectsize.jpeg",width = 800, height = 350)
ggplot(comps.MTSV,aes(estimate,trait))+geom_point(size=4,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)
dev.off()

jpeg("..//figure/USFS_effectsize.jpeg",width = 800, height = 350)
ggplot(comps.USFS,aes(estimate,trait))+geom_point(size=4,aes(color=category),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category))+geom_vline(aes(xintercept=0))+theme_base(base_size=11)+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)
dev.off()



###does the




###Analysis 4cd
source("..//Scripts/continuous_mod_prep.R")
HF<-filter(HF,species!=("QUAL"))
spforcontmods<-traits$species

HF.continuous.data<-filter(HF,species %in% c(spforcontmods))
HF.continuous.data<-left_join(HF.continuous.data,traits, by="species") ###This is the data for the continuous models

##zscore predictors for these models
HF.continuous.data$cent_pol<-(HF.continuous.data$pol-mean(HF.continuous.data$pol,na.rm=TRUE))/(2*sd(HF.continuous.data$pol,na.rm=TRUE))
HF.continuous.data$cent_minP<-(HF.continuous.data$min_precip-mean(HF.continuous.data$min_precip))/(2*sd(HF.continuous.data$min_precip))
HF.continuous.data$cent_floday<-(HF.continuous.data$fopn.jd-mean(HF.continuous.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.continuous.data$fopn.jd,na.rm=TRUE))

###phylo lm can only handle one value per species in the tree
pol.sum<- HF.continuous.data %>% group_by(name) %>%summarise(mean_pol=mean(cent_pol,na.rm=TRUE))
flo.sum<- HF.continuous.data %>% group_by(name)%>%summarise(mean_floday=mean(cent_floday,na.rm=TRUE))
P.sum<- HF.continuous.data %>% group_by(name)%>%summarise(mean_minP=mean(cent_minP,na.rm=TRUE))
funct.sum<-HF.continuous.data %>% group_by(name)%>%summarise(mean_funct_offset=mean(funct.offset,na.rm=TRUE))
phys.sum<-HF.continuous.data %>% group_by(name)%>%summarise(mean_phys_offset=mean(phys.offset,na.rm=TRUE))
mean.continuous.dat<-left_join(funct.sum,phys.sum)
mean.continuous.dat<-left_join(mean.continuous.dat,P.sum)
mean.continuous.dat<-left_join(mean.continuous.dat,flo.sum)
mean.continuous.dat<-left_join(mean.continuous.dat,pol.sum)


df.cont<-mean.continuous.dat[match(mytree.names, mean.continuous.dat$name),]
df.cont$funct.bin<-ifelse(df.cont$mean_funct_offset>=0,1,0)
df.cont$phys.bin<-ifelse(df.cont$mean_phys_offset>=0,1,0)
df.cont<-as.data.frame(df.cont)
###HF phylosig
d3<-comparative.data(HF.tree.pruned,df.cont,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
#physiological hysteranthy
phylo.d(d3, binvar=phys.bin) 
phylo.d(d3, binvar=funct.bin)

df.cont<- df.cont%>% remove_rownames %>% column_to_rownames(var="name")


write.csv(df.cont,"HarvardForest/HF.means.data.csv")
write.tree(HF.tree.pruned,"HarvardForest/HFtree.tre")
### continuous models
HF.cont.funct<-phylolm(mean_funct_offset~mean_pol+mean_floday+mean_minP+mean_pol:mean_floday+mean_pol:mean_minP+mean_floday:mean_minP,data=df.cont,phy=HF.tree.pruned,boot=599)
summary(HF.cont.funct)

HF.cont.phys<-phylolm(mean_phys_offset~mean_pol+mean_floday+mean_minP+mean_pol:mean_floday+mean_pol:mean_minP+mean_floday:mean_minP,data=df.cont,phy=HF.tree.pruned,boot=599)
summary(HF.cont.phys)

HF.fuct.continuous<-full_join(extract_coefs(HF.cont.funct),extract_CIs(HF.cont.funct),by="trait")
colnames(HF.fuct.continuous)<-c("trait","estimate","low","high")
HF.fuct.continuous$class<-"functional"

HF.phys.continuous<-full_join(extract_coefs(HF.cont.phys),extract_CIs(HF.cont.phys),by="trait")
colnames(HF.phys.continuous)<-c("trait","estimate","low","high")
HF.phys.continuous$class<-"physiological"


HF.cont.comps<-rbind(HF.phys.continuous,HF.fuct.continuous)
HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_pol")] <- "main effect: pollination syndrome"
HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_floday")] <- "main effect: flowering time"
HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_minP")] <- "main effect: minimum precipitation"
HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_floday:mean_minP")] <- "interaction: flowering x precip."
HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_pol:mean_minP")] <- "interaction: pollination x precip."
HF.cont.comps$trait[which(HF.cont.comps$trait=="mean_pol:mean_floday")] <- "interaction: pollination x flowering"

HF.cont.comps<-filter(HF.cont.comps,trait!="sigma2")

jpeg("..//figure/HF_cont_effectsize.jpeg",width = 800, height = 350)
ggplot(HF.cont.comps,aes(estimate,trait))+geom_point(size=4,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))
dev.off()

###binary doesn't really work in phylglm
HF.bin.funct<-phyloglm(funct.bin~mean_pol+mean_floday+mean_minP,df.cont, HF.tree.pruned, method = "logistic_MPLE", btol = 10, log.alpha.bound = 8,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)
summary(HF.bin.funct)

HF.bin.phys<-phyloglm(phys.bin~mean_pol+mean_floday+mean_minP,df.cont, HF.tree.pruned, method = "logistic_MPLE", btol = 800, log.alpha.bound = 10,
                      start.beta=NULL, start.alpha=NULL,
                      boot=599,full.matrix = TRUE)
summary(HF.bin.phys)

HF.fuct.bin<-full_join(extract_coefs(HF.bin.funct),extract_CIs(HF.bin.funct),by="trait")
colnames(HF.fuct.bin)<-c("trait","estimate","low","high")
HF.fuct.bin$class<-"functional"

HF.phys.bin<-full_join(extract_coefs(HF.bin.phys),extract_CIs(HF.bin.phys),by="trait")
colnames(HF.phys.bin)<-c("trait","estimate","low","high")
HF.phys.bin$class<-"physiological"


HF.bin.comps<-rbind(HF.phys.bin,HF.fuct.bin)
ggplot(HF.bin.comps,aes(estimate,trait))+geom_point(size=4,aes(color=class),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=class))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))


inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

###continuous models
modelcont.funct.noint <- brm(funct.offset~ cent_pol+cent_floday+cent_minP+cent_pol:cent_floday+cent_pol:cent_minP+cent_floday:cent_minP+(1|name), data = HF.continuous.data, 
                             family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelcont.phys.noint <- brm(phys.offset~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                            family = gaussian(), cov_ranef = list(name= A),iter=3000) 

modelbin.funct.noint <- brm(hyst.funct~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                            family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) ##1 divergent transition

modelbin.phys.noint <- brm(hyst.phys~ cent_pol+cent_floday+cent_minP+(1|name), data = HF.continuous.data, 
                           family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=3000) 

#extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.10,.25,.75,0.90))),"trait")
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
intra.df$offset.cent<-intra.df$offset-mean(intra.df$offset,na.rm=TRUE)

####now run each species model's seperately
df.intra.alnus<-filter(intra.df,taxa=="Alnus glutinosa")
df.intra.frax<-filter(intra.df,taxa=="Fraxinus excelsior")
df.intra.bet<-filter(intra.df,taxa=="Betula pendula")
df.intra.aes<-filter(intra.df,taxa=="Aesculus hippocastenum")




lm(offset.cent~soil.cent,data=df.intra.alnus)
lm(offset~soil.cent,data=df.intra.aes)
lm(offset~soil.cent,data=df.intra.frax)

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

library(lme4)
summary(lmer(offset~soil.cent+(1|taxa),data=intra.df))
summary(lmer(offset~flo.cent+soil.cent+(1|taxa),data=intra.df))

jpeg("..//figure/SM_comp.jpeg")
par(mfrow=c(1, 1) )
plot(c(-10,10), c(0,6), type = "n", xlab = "Soil moisture (deviation from mean)", ylab = "change in FLS offset", bty='l')
segments(x=0,y0=4.809891 ,x1=10,y1=4.809891-0.036298*10,col="black",lty="solid",lwd=3 )
segments(x=0,y0=4.809891 ,x1=-10,y1=4.809891+0.036298*10,col="black",lty="solid" ,lwd=3)
segments(x=0,y0=1.799010 ,x1=10,y1=1.799010+0.133845*10,col="black",lty="dashed",lwd=3 )
segments(x=0,y0=1.799010 ,x1=-10,y1=1.799010-0.133845*10,col="black",lty="dashed",lwd=3 )
dev.off()


priorz<-get_prior(offset~flo.cent+leaf.cent+soil.cent,data=df.intra.alnus)
mod.intra.flo.sm.aln<-brm(offset~flo.cent+soil.cent,data=df.intra.alnus,prior=priorz)
summary(mod.intra.flo.sm.aln)

priorz.fr<-get_prior(offset~flo.cent+soil.cent,data=df.intra.frax)
mod.intra.flo.sm.frax<-brm(offset~flo.cent*soil.cent,data=df.intra.frax,prior=priorz)
summary(mod.intra.flo.sm.frax)
mod.intra.flo.sm.sm<-brm(offset~flo.cent*soil.cent,data=df.intra.bet,prior=priorz)
summary(mod.intra.flo.sm.sm)

alnus.itra<-extract_coefs4HF(mod.intra.flo.sm.aln)
frax.itra<-extract_coefs4HF(mod.intra.flo.sm.frax)
bet.itra<-extract_coefs4HF(mod.intra.flo.sm.sm)
save.image("paper_full_analysis.RData")
