
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

load("consensis_plot.Rdata")

### For Final plot
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

######### phyloglm requires species names to be in rownames
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
silv.data<- silv.data %>% remove_rownames %>% column_to_rownames(var="name")

###flip the sign of floday so it can be interpreted better on graph
##MTSV and USFS with Interactions
mich.data$flo_cent.neg<--(mich.data$flo_cent)
silv.data$flo_cent.neg<--(silv.data$flo_cent)

z.funct.drought<-phyloglm(pro2~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=599,full.matrix = TRUE)


z.phys.drought<-phyloglm(pro3~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg,mich.data, mich.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                         start.beta=NULL, start.alpha=NULL,
                         boot=599,full.matrix = TRUE)

z.funct.drought.silvics<-phyloglm(pro2~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                  start.beta=NULL, start.alpha=NULL,
                                  boot=599,full.matrix = TRUE)

z.phys.drought.silvics<-phyloglm(pro3~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg,silv.data, silv.tre, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                                 start.beta=NULL, start.alpha=NULL,
                                 boot=599,full.matrix = TRUE)

###combine MTSV and USFS model out put
##prep for plotting
mich.funct.dat<-full_join(extract_coefs(z.funct.drought),extract_CIs(z.funct.drought),by="trait")
colnames(mich.funct.dat)<-c("trait","estimate","low","high")
mich.funct.dat$class<-"functional-MTSV"

mich.phys.dat<-full_join(extract_coefs(z.phys.drought),extract_CIs(z.phys.drought),by="trait")
colnames(mich.phys.dat)<-c("trait","estimate","low","high")
mich.phys.dat$class<-"physiological-MTSV"

silv.funct.dat<-full_join(extract_coefs(z.funct.drought.silvics),extract_CIs(z.funct.drought.silvics),by="trait")
colnames(silv.funct.dat)<-c("trait","estimate","low","high")
silv.funct.dat$class<-"functional-USFS"

silv.phys.dat<-full_join(extract_coefs(z.phys.drought.silvics),extract_CIs(z.phys.drought.silvics),by="trait")
colnames(silv.phys.dat)<-c("trait","estimate","low","high")
silv.phys.dat$class<-"physiological-USFS"

MTSV<-rbind(mich.phys.dat,mich.funct.dat)
USFS<-rbind(silv.phys.dat,silv.funct.dat)

comps<-rbind(MTSV,USFS)
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

##tempory plot of these two only
pd=position_dodgev(height=0.4)
ggplot(comps,aes(estimate,trait))+geom_point(size=4,aes(color=category,shape=data),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=low,xmax=high,color=category,linetype=data))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))+xlim(-8,8)

####add HF
HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
HF$phys.offset<-HF$bb.jd-HF$fbb.jd
HF$funct.offset<-HF$l75.jd-HF$fopn.jd

source("..//Scripts/continuous_mod_prep.R")
HF<-filter(HF,species!=("QUAL"))
spforcontmods<-df$species

HF.continuous.data<-filter(HF,species %in% c(spforcontmods))
HF.continuous.data<-left_join(HF.continuous.data,traits, by="species") ###This is the data for the continuous models
##zscore predictors for these models
HF.continuous.data$pol_cent<-(HF.continuous.data$pol-mean(HF.continuous.data$pol,na.rm=TRUE))/(2*sd(HF.continuous.data$pol,na.rm=TRUE))
HF.continuous.data$precip_cent<-(HF.continuous.data$min_precip-mean(HF.continuous.data$min_precip))/(2*sd(HF.continuous.data$min_precip))
HF.continuous.data$flo_cent<-(HF.continuous.data$fopn.jd-mean(HF.continuous.data$fopn.jd,na.rm=TRUE))/(2*sd(HF.continuous.data$fopn.jd,na.rm=TRUE))
HF.continuous.data$flo_cent.neg<--(HF.continuous.data$flo_cent)

inv.phylo <- MCMCglmm::inverseA(HF.tree.pruned, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

modelcont.funct <- brm(funct.offset~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.continuous.data, 
                             family = gaussian(), cov_ranef = list(name= A),iter=4000) 


modelcont.phys <- brm(phys.offset~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.continuous.data, 
                            family = gaussian(), cov_ranef = list(name= A),iter=4000) 

modelbin.funct<- brm(hyst.funct~ pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.continuous.data, 
                            family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000) 

modelbin.phys<- brm(hyst.phys~pol_cent+flo_cent.neg+precip_cent+precip_cent:flo_cent.neg+precip_cent:pol_cent+pol_cent:flo_cent.neg+(1|name), data = HF.continuous.data, 
                           family = bernoulli(link="logit"), cov_ranef = list(name= A),iter=4000) 

#pp_check(modelcont.funct)
#pp_check(modelcont.phys)
#pp_check(modelbin.funct)
#pp_check(modelbin.phys)

extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.975))),"trait")
}

funct.cont<-extract_coefs4HF(modelcont.funct)
funct.cont$class<-"functional"
funct.bin<-extract_coefs4HF(modelbin.funct)
funct.bin$class<-"functional"

phys.cont<-extract_coefs4HF(modelcont.phys)
phys.bin<-extract_coefs4HF(modelbin.phys)

phys.cont$class<-"physiological"
phys.bin$class<-"physiological"

cont<-rbind(phys.cont,funct.cont)
bin<-rbind(phys.bin,funct.bin)
cont$data_type<-"continuous"
bin$data_type<-"binary"

hfboth<-rbind(cont,bin)    

pd=position_dodgev(height=0.4)
ggplot(hfboth,aes(Estimate,trait))+geom_point(size=4,aes(color=class,shape=(data_type)),position=pd)+geom_errorbarh(position=pd,width=0,aes(xmin=Q2.5,xmax=Q97.5,color=class,linetype=data_type))+geom_vline(aes(xintercept=0))+theme_base(base_size = 11)+scale_color_manual(values=c("orchid4", "springgreen4"))


colnames(hfboth)[colnames(hfboth)=="class"] <- "category"
colnames(hfboth)[colnames(hfboth)=="Estimate"] <- "estimate"

comps2<-comps
colnames(comps2)<-c("trait","estimate","Q2.5","Q97.5","class","category","data")
comps2<-dplyr::select(comps2,-class)
#hfboth<-dplyr::select(hfboth,-Q25)
#hfboth<-dplyr::select(hfboth,-Q75)
hfboth<-dplyr::select(hfboth,-Est.Error)

colnames(hfboth)
colnames(comps2)
hfboth$data<-"HF"
comps2$data_type<-"binary"

comps2<-comps2[,c(1,2,3,4,5,7,6)]

alleffectos<-rbind(comps2,hfboth)
alleffectos<-dplyr::filter(alleffectos,trait!="(Intercept)")
alleffectos<-dplyr::filter(alleffectos,trait!="Intercept")
alleffectos$trait[which(alleffectos$trait=="pol_cent")]<-"pollination syndrome"
alleffectos$trait[which(alleffectos$trait=="flo_cent.neg")]<- "earlier flowering"
alleffectos$trait[which(alleffectos$trait=="precip_cent")]  <- "water dynamics"
alleffectos$trait[which(alleffectos$trait=="pol_cent:precip_cent")]<- "pollination:water dynamics"
alleffectos$trait[which(alleffectos$trait=="pol_cent:flo_cent.neg")]<-"pollination:flowering"
alleffectos$trait[which(alleffectos$trait=="flo_cent.neg:precip_cent")]<-"flowering:water dynamics"




allos<-filter(alleffectos,data!="HF")

jpeg("..//figure/option1.jpeg",width = 8.6, height = 4, units = 'in', res=200)
pd=position_dodgev(height=0.4)
alleffectos %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=data_type,color=data),position=pd,size=3,stroke=1.5)+scale_shape_manual(name="data type",values=c(15,16))+scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("dotted","dotted"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-30,30)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+guides(size = "legend", linetype= "none")+
  annotate("text", x = 24.2, y = 6.4, label = "Hysteranthy",fontface="bold",size=3)+annotate("text", x = -24.9, y = 6.4, label = "Seranthy",fontface="bold",size=3)+
  #guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))
  facet_wrap(~category)
dev.off()

binneys<-filter(alleffectos, data_type=="binary")

pd=position_dodgev(height=0.4)
twoa<-binneys %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(color=data,fill=category),shape=21,position=pd,size=3,stroke=1.5)+
 scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=data,linetype=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("solid","solid"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-10,12)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+guides(size = "legend", linetype= "none")+
  annotate("text", x = 11, y = 6.4, label = "Hysteranthy",fontface="bold",size=3)+annotate("text", x = -9, y = 6.4, label = "Seranthy",fontface="bold",size=3)+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))


cont<-filter(alleffectos, data_type!="binary")

pd=position_dodgev(height=0.4)
twob<-cont %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(color=data,fill=category),shape=22,position=pd,size=3,stroke=1.5)+
  scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=data,linetype=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("solid","solid"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-30,30)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+guides(size = "legend", linetype= "none")+
  annotate("text", x = 21.2, y = 6.4, label = "Hysteranthy",fontface="bold",size=3)+annotate("text", x = -21.9, y = 6.4, label = "Seranthy",fontface="bold",size=3)+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))

jpeg("..//figure/option2.jpeg",width = 12, height = 4, units = 'in', res=100)
ggpubr::ggarrange(twoa,twob)
dev.off()

pd=position_dodgev(height=0.4)
threea<-allos %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=data_type,color=data,fill=category),position=pd,size=3)+scale_shape_manual(name="data type",values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data,alpha=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(name=NULL, values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-8,8)+
  scale_color_manual(name="data set",values=c("darkgoldenrod1", "springgreen4"))+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))+
  guides(size = "legend", linetype= "none")+
  guides(size = "legend", shape= "none")+
  annotate("text", x = 6.4, y = 6.5, label = "Hysteranthy",fontface="bold")+annotate("text", x = -7, y = 6.5, label = "Seranthy",fontface="bold")

HFer<-filter(alleffectos,data=="HF")
threeb<-HFer %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=data_type,color=data,fill=category),position=pd,size=3,stroke=1.5)+scale_shape_manual(name="data type",values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data,fill=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("solid","solid"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-30,30)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+guides(size = "legend", linetype= "none")+
  annotate("text", x = 24.2, y = 6.4, label = "Hysteranthy",fontface="bold",size=3)+annotate("text", x = -24.9, y = 6.4, label = "Seranthy",fontface="bold",size=3)+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))+
 guides(size = "legend", color= "none")

jpeg("..//figure/option3.jpeg",width = 12, height = 4, units = 'in', res=100)
ggpubr::ggarrange(threea,threeb)
dev.off()



foura<-binneys %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(shape=data_type,color=data,fill=category),position=pd,size=3,stroke=1.5)+scale_shape_manual(name="data type",values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data,alpha=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(name=NULL, values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-8,12)+
  facet_wrap(~data)+scale_color_manual(name="data set",values=c("orchid4","darkgoldenrod1", "springgreen4"))+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))+
  guides(size = "legend", linetype= "none")+
  annotate("text", x = 26, y = 6.5, label = "Hysteranthy",fontface="bold")+annotate("text", x = -27.5, y = 6.5, label = "Seranthy",fontface="bold")
jpeg("..//figure/option4.jpeg",width = 12, height = 6, units = 'in', res=100)
ggpubr::ggarrange(foura,twob,nrow=2)
dev.off()
alleffectos.noint<-dplyr::filter(alleffectos,trait %in% c("pollination syndrome"    ,   "earlier flowering"     ,     "water dynamics"))

jpeg("..//figure/option5.jpeg",width = 8.6, height = 4, units = 'in', res=200)
pd=position_dodgev(height=0.6)
alleffectos.noint %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
  ggplot(aes(estimate,trait))+geom_point(aes(color=data,shape=data_type,fill=category),position=pd,size=3,stroke=1.5)+
  scale_shape_manual(name="data type",values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=data,linetype=data_type,fill=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(values=c("solid","solid"))+theme_base(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-15,30)+
  scale_color_manual(values=c("orchid4","darkgoldenrod1", "springgreen4"))+guides(size = "legend", linetype= "none")+
  annotate("text", x = 24.2, y = 3.4, label = "Hysteranthy",fontface="bold",size=3)+annotate("text", x = -9.9, y = 3.4, label = "Seranthy",fontface="bold",size=3)+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",physiological="gray"))))

dev.off()










####now add intra specific
aln10<-read.csv("datasheets_derived/PEP/alnus10_delta_hyst.csv",header=TRUE)
frax10<-read.csv("datasheets_derived/PEP/fraxinus10_delta_hyst.csv",header=TRUE)
bet10<-read.csv("datasheets_derived/PEP/betpen10_delta_hyst.csv",head=TRUE)
aes10<-read.csv("datasheets_derived/PEP/aes10_delta_hyst.csv",header=TRUE)

aln10$taxa<-"Alnus glutinosa" ##assign species
frax10$taxa<-"Fraxinus excelsior"
aes10$taxa<-"Aesculus hippocastenum"

intra.d<-rbind(aln10,frax10,bet10,aes10)
moist.aug<-raster("climate_data/grids_germany_multi_annual_soil_moist_1991-2010_08.asc") ##Gause Kruger 3 for August
moist.apr<-raster("climate_data/grids_germany_multi_annual_soil_moist_1991-2010_04.asc")

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
intra.df$flo.cent<-intra.df$flower-mean(intra.df$flower,na.rm=TRUE)/(2*sd(intra.df$flower,na.rm=TRUE))
intra.df$soil.cent<-intra.df$SM-mean(intra.df$SM,na.rm=TRUE)/(2*sd(intra.df$SM,na.rm=TRUE))
intra.df$flo.cent.neg<--(intra.df$flo.cent)

##try it without aesculus
intra.df.hystonly<-dplyr::filter(intra.df, taxa!="Aesculus hippocastenum")
smpriorz<-brms::get_prior(offset~flo.cent.neg*soil.cent,data=intra.df.hystonly) #The interaction term here was 0


smflo.bayes<-brm(offset~flo.cent.neg*soil.cent,data=intra.df.hystonly,prior=smpriorz)
summary(smflo.bayes)
pepbayes<-extract_coefs4HF(smflo.bayes)
colnames(alleffectos)
colnames(pepbayes)
pepbayes<-dplyr::select(pepbayes,-Est.Error)
pepbayes$category<-"other"
pepbayes$data_type<-"continuous"
pepbayes$data<-"PEP725"
colnames(pepbayes)<-c("trait"  ,   "estimate",  "Q2.5"    ,  "Q97.5"   ,  "category",  "data_type" ,"data"  )
pepbayes<-dplyr::filter(pepbayes,trait!="Intercept")
head(pepbayes)
head(alleffectos)
pepbayes$trait[which(pepbayes$trait=="flo.cent.neg")] <- "earlier flowering"
pepbayes$trait[which(pepbayes$trait=="soil.cent")] <- "water dynamics"
pepbayes$trait[which(pepbayes$trait=="flo.cent.neg:soil.cent")] <- "flowering:water dynamics"

alleffectos$trait[which(alleffectos$trait=="pol_cent")]<-"pollination syndrome"
alleffectos$trait[which(alleffectos$trait=="flo_cent.neg")]<- "earlier flowering"
alleffectos$trait[which(alleffectos$trait=="precip_cent")]  <- "water dynamics"
alleffectos$trait[which(alleffectos$trait=="pol_cent:precip_cent")]<- "pollination:water dynamics"
alleffectos$trait[which(alleffectos$trait=="pol_cent:flo_cent.neg")]<-"pollination:flowering"
alleffectos$trait[which(alleffectos$trait=="flo_cent.neg:precip_cent")]<-"flowering:water dynamics"

biggest<-rbind(alleffectos,pepbayes)
pd=position_dodgev(height=0.4)



jpeg("..//figure/allmods_effectsizes_combined.jpeg",width = 8, height = 6, units = 'in', res=300)
biggest %>%
  arrange(estimate) %>%
  mutate(trait = factor(trait, levels=c("flowering:water dynamics","pollination:flowering","pollination:water dynamics","earlier flowering","water dynamics","pollination syndrome"))) %>%
ggplot(aes(estimate,trait))+geom_point(aes(shape=data_type,color=data,fill=category),position=pd,size=3,stroke=1.5)+scale_shape_manual(name="data type",values=c(21,22))+scale_fill_manual(values=c(functional="black",physiological="grey",other="beige"))+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,linetype=data_type,color=data,alpha=category),position=pd,width=0)+scale_alpha_manual(values=c(1,1,1))+
  scale_linetype_manual(name=NULL, values=c("solid","solid"))+theme_linedraw(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+xlim(-30,30)+
  scale_color_manual(name="data set",values=c("orchid4","darkgoldenrod1", "springgreen4","skyblue1"))+
  guides(fill=guide_legend(override.aes=list(colour=c(functional="black",other="beige",physiological="gray"))))+
  guides(size = "legend", linetype= "none")+
  annotate("text", x = 26, y = 6.5, label = "Hysteranthy",fontface="bold")+annotate("text", x = -27.5, y = 6.5, label = "Seranthy",fontface="bold")

dev.off()

library("ggeffects")
##take a model where there are interactions

goober<- ggpredict(modelbin.phys,c("precip_cent","pol_cent","flo_cent.neg[.2]"), ci.lvl=0.50) ##May 15
apc.phys<-plot(goober)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1),labels=c(15,33,50,67,84,101))+
  xlab("Min. precipitation across range (cm)")+ylab("Likelihood of flower before leaf budburst")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag="a)")+theme_linedraw()

goober2<- ggpredict(modelcont.funct,c("precip_cent","pol_cent","flo_cent.neg[.2]"), ci.lvl=0.50)  #may 15
apc.funct<-plot(goober2)+scale_x_continuous(breaks =c(-1.5,-1.0,-0.5,0,0.5,1),labels=c(15,33,50,67,84,101))+
  xlab("Min. precipitation across range (cm)")+ylab("Flowers opening to 75% leaf expansion (days)")+scale_colour_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+scale_fill_manual(name="pollination syndrome",labels=c("biotic","wind"),values=c("coral4","royalblue2"))+
  labs(title = NULL,tag="b)")+theme_linedraw()

jpeg("..//figure/HF_apcs.jpeg",width = 8, height = 6, units = 'in', res=300)
ggpubr::ggarrange(apc.phys,apc.funct,ncol=2,legend="bottom",common.legend = TRUE)
dev.off()
### for insect offset decreases as minp goes up, but for wind chance of hysteranthy goes up
z.p<-c(-1.5,-1.0,-0.5,0,0.5,1)

(1*(2*sd(HF.continuous.data$min_precip,na.rm=TRUE))+mean(HF.continuous.data$min_precip,na.rm=TRUE))*2.54

.2*(2*sd(HF.continuous.data$fbb.jd,na.rm=TRUE))+mean(HF.continuous.data$fbb.jd,na.rm=TRUE)
31+28+31+31+15
save.image("consensis_plot.Rdata")

citation("base")
"citation"

##### do a quicky MTSV and USFS Comparison

