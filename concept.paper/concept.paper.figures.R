setwd("~/Documents/git/proterant/concept.paper/")
library(ggthemes)
source("concept.paper.analysis.R")
###run concept.paper.analysis.R all the way line by line. Sources yeilds an error
#Figure 1: Harvard forest hysteranthy changes
###functional hysteranthy for harvard forest
source("HF.hyst.4.manuscript.R")
fun
phys
####TO DO
#make species names more visable
#customize symbol
#elimiate entries missing a phenophase?
###figure 2: tree with hysteranthy information
source("make_circle_tree.R")
p


###Figure 3:effect plots
##3a full dataset
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

functplot<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme_tufte()+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
functplot
#3b restricted early data
source("earlyfloweronly.R") ### I think you have to run this line by line, sourcing drops all the tips
bootest<-as.data.frame(earlycent$coefficients)
bootconf<-as.data.frame(earlycent$bootconfint95)
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

earlplot<-ggplot(bootmich,aes(estimate,trait))+geom_point(size=2.5)+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+theme_tufte()+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+theme(axis.text = element_text(size=14, hjust = .5))+guides(color="none")
earlplot

###Suppliment
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

physplot<-ggplot(bootmich,aes(estimate,trait))+geom_point()+geom_segment(aes(y=trait,yend=trait,x=low,xend=high))+ggtitle("Physiological Hysteranthy")+theme_tufte()+theme(panel.border=element_rect(aes(color=blue)))+geom_vline(aes(xintercept=0,color="red"))+theme(plot.title = element_text(hjust = 0.5))+guides(color="none")
physplot
###
source("make_hypothetical.R")
citation("phylolm")
citation("caper")
citation
?phylo.d()
