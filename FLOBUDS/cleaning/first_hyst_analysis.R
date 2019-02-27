#####First flobuds analysis. Question: How does first flowering (60), change relative to first leafout day (15)?
# Begun by Dan on April 13 2018
setwd("~/Documents/git/proterant/FLOBUDS")
source("cleaning/cleaning.R")

#############CREATE A DATA SHEET THAT HAS THE THREE FLOWER TYPES IN SEPARATE COLUMNS#########################
#################THIS WILL BE USEFUL FOR COMPARING CHANGES IN PROTANDRY OR PROTOGYNY######################

####find the first day when species reached 15
d.leaf<-filter(dx,leafphase==15)
first<-aggregate(d.leaf$doy.final, by = list(d.leaf$id), min)
####combine with all data
dater<-as.data.frame(unique(dx$id))
colnames(dater)<- c("id")
colnames(first)<-c("id","leaf_day(15)")
dater<-full_join(dater,first,by="id") ###now you have a data set with first leaves



### flowers (currrently mixed only)
unique(dx$flophase)
d.flo<-filter(dx,flophase %in% c("60","60-F","60-M")) ###This is all the 60s ever I have
firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),min )

colnames(firstflo)<-c("id","flo_day")
dater<-full_join(dater,firstflo,by="id")
###add female
#d.floF<-filter(dx,flophase=="60-F")
#firstfloF<-aggregate(d.floF$doy.final, by = list(d.floF$id), min) ##This stage is problematc

#colnames(firstfloF)<-c("id","flo_dayF")
#dater<-full_join(dater,firstfloF,by="id")
###add male
#d.floM<-filter(dx,flophase=="60-M")
#firstfloM<-aggregate(d.floM$doy.final, by = list(d.floM$id), min)
#colnames(firstfloM)<-c("id","flo_dayM")
#dater<-full_join(dater,firstfloM,by="id")
##### pull experimenta; information to merge. THis give a full data sheet 
###in which male, female and mixd flowering can be compared
dxx<-dplyr::select(dx,1:7)
dxx<-distinct(dxx)
good.dat<-left_join(dater,dxx)
###make a dichogamy sheet
#mono<-filter(good.dat, GEN.SPA %in% c( "COM.PER","COR.COR"))
#write.csv(mono,"dicogamy.csv")
###continue
good.dat<-gather(good.dat,sex,DOY,2:5)



##########################################################################################

############MAKE A DATE SHEET THAT AGGREGATES ALL FLOWER TO EVALUATE COARSE HYSTERANTHY############
unique(dx$flophase)
DAT<-separate(dx,flophase,c("abosolute_flower","sex_ind"),sep="-")
L<-filter(DAT,leafphase==9)
L1<-aggregate(L$doy.final, by = list(L$id), min)

datUM<-as.data.frame(unique(DAT$id))
colnames(datUM)<- c("id")
colnames(L1)<-c("id","Lbb_day(9)")
datUM<-full_join(datUM,L1,by="id") ###now you have a data set with first leaves

### flowers (currrently mixed only)
unique(DAT$abosolute_flower)
Fl<-filter(DAT,abosolute_flower==60)
Fl1<-aggregate(Fl$doy.final, by = list(Fl$id), min)
colnames(Fl1)<-c("id","flo_day(60)")
datUM<-full_join(datUM,Fl1,by="id")
great.dat<-left_join(datUM,dxx)

##add leeaf out
LLL<-filter(DAT,leafphase==15)
L3<-aggregate(LLL$doy.final, by = list(LLL$id), min)

datUM3<-as.data.frame(unique(DAT$id))
colnames(datUM3)<- c("id")
colnames(L3)<-c("id","leaf_day(15)")
great.dat<-left_join(great.dat,L3)
###add expansion
LL<-filter(DAT,leafphase==11)
L2<-aggregate(LL$doy.final, by = list(LL$id), min)

datUM2<-as.data.frame(unique(DAT$id))
colnames(datUM2)<- c("id")
colnames(L2)<-c("id","Lexpand_day(11)")
great.dat<-left_join(great.dat,L2)

BB<-filter(DAT,leafphase==07)
LBB<-aggregate(BB$doy.final, by = list(BB$id), min)

datUM2<-as.data.frame(unique(DAT$id))
colnames(datUM2)<- c("id")
colnames(L2)<-c("id","Lexpand_day(11)")
great.dat<-left_join(great.dat,L2)


ncol(great.dat)
great.dat<-great.dat[,c(1,4,5,6,7,8,9,3,2,10,11)]


###add survival analysis columns
surv<-read.csv("data/flobuds.eval.csv",header=TRUE)
colnames(surv)
colnames(great.dat)
surv<-unite(surv,treatcode,Force,Light,Chill,sep="_",remove = FALSE)
surv$Light<-ifelse(surv$treatcode=="W_S_0","L",surv$Light)
surv$Light<-ifelse(surv$treatcode=="W_L_0","S",surv$Light)
###clean this

head(great.dat)
surv<-unite(surv,id,name,good_flaskid,sep="_")
intersect(surv$id,great.dat$id) 
datty<-left_join(great.dat,surv)
##a data sheet for first analysis
write.csv(datty,"first.event.dat.csv") 
#great.dat<-unite(great.dat,treatment,Force,Light,Chill,sep="" )
#####################################
###good.dat, great.dat are the files to analyze
########################################################################
### metrics for great date##############################################
###how many entries have full entries?
full<-subset(great.dat, !is.na(great.dat$bb_day)&!is.na(great.dat$flo_day))
nrow(full)#129 our of 576
table(full$GEN.SPA)
table(full$treatment)
### how many have flowering?
flowerfun<-subset(great.dat,!is.na(great.dat$flo_day))
nrow(flowerfun) #142 out of 576
table(flowerfun$GEN.SPA)
table(flowerfun$treatment)

#one or the other:
something<-subset(great.dat, !is.na(great.dat$leaf_day)|!is.na(great.dat$flo_day))
som<-gather(something,phenophase,DOY,2:3)

nrow(something) #353
353/576
table(something$GEN.SPA)
table(something$treatment)

############################################################################################################
###plot the raw data

an.data<-gather(great.dat,Phenophase,DOY,2:3)
unique(an.data$GEN.SPA)
##exclude betula, ace sac, and
bigsp<-filter(an.data, GEN.SPA %in% c("ACE.PEN","ACE.RUB","PRU.VIR", "COM.PER","COR.COR","ILE.MUC","ILE.VER","VIB.ACE", "PRU.PEN","VAC.COR"))
ggplot(bigsp, aes(x=treatment, y=DOY,color=Phenophase))+stat_summary()+geom_point(size=.25)+facet_wrap(~GEN.SPA)


###raw data-How much data do we have 
x<-subset(great.dat, !is.na(flo_day))
y<-subset(great.dat,!is.na(leaf_day))
table(x$GEN.SPA)
table(y$GEN.SPA)
xx<-rbind(x,y)
xx<-gather(xx,phenophase,DOY,2:3)
p<-ggplot(x,aes(treatment))+stat_count(color="pink",geom="bar")+facet_wrap(~GEN.SPA)
pp<- ggplot(y,aes(treatment))+stat_count(color="green",geom="bar")+facet_wrap(~GEN.SPA)
ggplot(p+pp)
library(gridExtra)
grid.arrange(p, pp)

####simple plots

ggplot(great.dat,aes(as.factor(Chill),flo_day))+geom_point(aes(color=Force,shape=Light))+facet_wrap(~GEN.SPA)+ggtitle("First Flowering")
ggplot(great.dat,aes(as.factor(Chill),leaf_day))+geom_point(aes(color=Force,shape=Light))+facet_wrap(~GEN.SPA)+ggtitle("First Leafout")


###basic linear model

ggplot(great.dat,aes(flo_day))+geom_density() ### not so normal
ggplot(great.dat,aes(leaf_day))+geom_density() ### not so normal


