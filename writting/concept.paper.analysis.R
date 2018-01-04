####concept paper code: THis is based off hyst_final_analysis.R. but only what we need to the manuscript
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")
library("ape")
library("phytools")
library("geiger")
library("gbm")
library("pez")
library(caper)
library(picante)
library("tidyverse")
library(boot)
library("phylolm")
library("ggplot2")
library(arm)

#########READ IN ALL DATA AND ASSOCIATED TREES##################

mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")
###cleaning:
#clean av fruit time
mich.data$av_fruit_time[mich.data$av_fruit_time=="persistant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="persitant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="unreported"]<-9    
mich.data$av_fruit_time<-as.numeric(mich.data$av_fruit_time)

mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
#mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9                                      
mich.data$fruiting<-as.numeric(mich.data$fruiting)

mich.data["pro3"]<-NA
mich.data$pro3[mich.data$Phen.sequence == "pro"] <- 1
mich.data$pro3[mich.data$Phen.sequence == "pro/syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn/ser"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "ser"] <- 0 
mich.data$pro3[mich.data$Phen.sequence== "hyst"] <- 0

###phylo.D
set.seed(122)
mich.tree$node.label<-NULL
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloD <- phylo.d(d, binvar=pro) ###regular hysteranthy
PhyloD
##functionalhysteranthy
PhyloPro2<-phylo.d(d,binvar=pro2)
PhyloPro2
d<-comparative.data(mich.tree,mich.data,name,vcv = TRUE,vcv.dim = 2, na.omit = FALSE)
PhyloPro3<-phylo.d(d,binvar=pro3)
PhyloPro3
###centering
mich.data$height_cent<-(mich.data$heigh_height-mean(mich.data$heigh_height))/(2*sd(mich.data$heigh_height))
mich.data$fruit_cent<-(mich.data$fruiting-mean(mich.data$fruiting))/(2*sd(mich.data$fruiting))
mich.data$flo_cent<-(mich.data$flo_time-mean(mich.data$flo_time))/(2*sd(mich.data$flo_time))
mich.data$pol_cent<-(mich.data$pol-mean(mich.data$pol))/(2*sd(mich.data$pol))
mich.data$av_fruit_time_cent<-(mich.data$av_fruit_time-mean(mich.data$av_fruit_time))/(2*sd(mich.data$av_fruit_time))

###prepare for modeling
mich.data<-  mich.data %>% remove_rownames %>% column_to_rownames(var="name")
earl<-  earl%>% remove_rownames %>% column_to_rownames(var="name")
mich.data$height10<-mich.data$heigh_height/10
###models:
mich5<-phyloglm(pro2~pol+height10+flo_time+fruiting+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                start.beta=NULL, start.alpha=NULL,
                boot=50,full.matrix = TRUE)
summary(mich5)

Mich5cent.funct<-phyloglm(pro2~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)

summary(Mich5cent.funct)

Mich5cent.super<-phyloglm(pro3~pol+height_cent+flo_cent+fruit_cent+shade_bin,mich.data, mich.tree, method = "logistic_MPLE", btol = 100, log.alpha.bound = 10,
                          start.beta=NULL, start.alpha=NULL,
                          boot=50,full.matrix = TRUE)
summary(Mich5cent.super)

###effect plots
bootest<-as.data.frame(Mich5cent.funct$coefficients)
bootconf<-as.data.frame(Mich5cent.funct$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="alpha")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
bootmich$trait[bootmich$trait=="shade_bin"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="fruit_cent"]<-"fruit timing"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

functplot<-ggplot(bootmich,aes(estimate,trait))+geom_point()+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+ggtitle("Functional Hysteranthy")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+theme(plot.title = element_text(hjust = 0.5))+guides(color="none")

####now physiological
bootest<-as.data.frame(Mich5cent.super$coefficients)
bootconf<-as.data.frame(Mich5cent.super$bootconfint95)
bootconf<-as.data.frame(t(bootconf))

bootest<-rownames_to_column(bootest, "trait")
bootconf<-rownames_to_column(bootconf, "trait")
bootmich<-full_join(bootconf,bootest, by="trait")
colnames(bootmich)<-c("trait","low","high","estimate")
bootmich<-dplyr::filter(bootmich, trait!="alpha")
bootmich<-dplyr::filter(bootmich, trait!="(Intercept)")
###names
bootmich$trait[bootmich$trait=="shade_bin"]<-"shade tolerance"
bootmich$trait[bootmich$trait=="pol"]<-"pollination syndrome"
bootmich$trait[bootmich$trait=="height_cent"]<-"max height"
bootmich$trait[bootmich$trait=="fruit_cent"]<-"fruit timing"
bootmich$trait[bootmich$trait=="flo_cent"]<-"flower timing"

physplot<-ggplot(bootmich,aes(estimate,trait))+geom_point()+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+ggtitle("Physiological Hysteranthy")+theme_light()+geom_vline(aes(xintercept=0,color="red"))+theme(plot.title = element_text(hjust = 0.5))+guides(color="none")




###average predictive comparsons using uncentered data except for height because it wouldnt converge
beta<-coef(mich5)
hi<-1
lo<-0
###for [pollination syndrome]
beta[2]
delta<-invlogit(beta[1]+beta[2]*hi+beta[3]*mich.data$height10+beta[4]*mich.data$flo_time+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*lo+beta[3]*mich.data$height10+beta[4]*mich.data$flo_time+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)

print(mean(delta))

###for flowering
earl<-4
mid<-5
delta<-invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height10+beta[4]*earl+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)-
  invlogit(beta[1]+beta[2]*mich.data$pol+beta[3]*mich.data$height10+beta[4]*mid+beta[5]*mich.data$fruiting+beta[6]*mich.data$shade_bin)

print(mean(delta))




###inter anual variation example not looking great

ts<-read.csv("individual_phenometrics_data.csv",header=TRUE)
ts<-dplyr::select(ts,Individual_ID,Phenophase_ID,Common_Name,Phenophase_Description,First_Yes_Year,First_Yes_DOY)
ts$First_Yes_Year<-as.character(ts$First_Yes_Year)
ts<-filter(ts, First_Yes_Year!="2015")
unique(ts$Phenophase_Description)
ts.func<-filter(ts, Phenophase_Description %in% c("Flowers or flower buds","Increasing leaf size"))
table(ts.func$Phenophase_Description)
q<-ggplot(ts.func, aes(x=First_Yes_Year, y=First_Yes_DOY, color=Phenophase_Description)) +
  stat_summary()+labs(title="tree spotters", x="Year", y="Days since initiation")
q+facet_wrap(~Common_Name)

#Harvard forest
hf<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
hf<-filter(hf,year>2004)
hf$year<-as.integer(hf$year)
###plot the budburst data###
new<-rename(hf,lbb.jd=bb.jd)
head(new)
new<-dplyr::select(new,-l75.jd)
head(new)
new<-gather(new,phenophase,eventday,lbb.jd:fbb.jd)
###Filter out species with insufficient flowering records
new<-filter(new, species %in% c( "ACRU","BEAL","QURU"))
q<-ggplot(new, aes(x=year, y=eventday, color=phenophase)) +
  stat_summary()+labs(title="Flower and Leaf Budburst at Harvard Forest", x="Year", y="Days since initiation")
q+facet_wrap(~species)
##Now plot flowers open vs. l75
new2<-dplyr::select(hf,-fbb.jd)
head(new2)
new2<-gather(new2,phenophase,eventday,l75.jd:fopn.jd)
###Filter out species with insufficient flowering records
new2<-filter(new2, species %in% c("ACRU","BEAL","QURU"))
q<-ggplot(new2, aes(x=year, y=eventday, color=phenophase)) +
  stat_summary()+labs(title="Flowering and Leaf Expansion at Harvard Forest", x="Year", y="Days since initiation")
q+facet_wrap(~species)
