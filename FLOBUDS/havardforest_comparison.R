rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")
library(ggplot2)
library(tidyr)
library(dplyr)
library(chillR)
library(brms)
library(tibble)
library(lme4)
library("lmerTest")
library(RColorBrewer)
library(ggstance)

HF<-read.csv("HarvardForest/hf003-05-mean-ind.csv",header=TRUE)
#HF2<-read.csv("HarvardForest/hf003-06-mean-spp.csv",header=TRUE)



##Hypothesis 1: variation in FLS is a product of variation in climate between phases
#https://cran.r-project.org/web/packages/chillR/vignettes/hourly_temperatures.html
weather<-read.csv("..//FLOBUDS/data/hf000-01-daily-m.csv",header = TRUE)
weather<-dplyr::select(weather,c("date","airtmax","airtmin"))
weather<-separate(weather,date,c("Year","Month","Day"),sep="-",remove=TRUE)
colnames(weather)<-c("Year","Month","Day","Tmax","Tmin")
sapply(weather,mode) #mode(weather)
weather$Year<-as.numeric(weather$Year)
weather$Month<-as.numeric(weather$Month)
weather$Day<-as.numeric(weather$Day)

unique(weather$Year)
weather<-filter(weather,Year>=1989)

#all_daylengths<-cbind(JDay=1:365,sapply(daylength(latitude=42.5,JDay=1:365),cbind)) ## calculate day length at HF on every day of year
#ad<-as.data.frame(all_daylengths) ## data frame of day length


weather<-make_all_day_table(weather)

hourtemps<-stack_hourly_temps(weather, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)


HF.bb<-filter(HF,!is.na(bb.jd))
indies<-unique(HF.bb$tree.id)
df<-data.frame(Season=character(),End_year=numeric(),Season_days=numeric(),Data_days=numeric(), Perc_complete=numeric(),Chilling_Hours=numeric(),Utah_Model=numeric(),Chill_portions=numeric(), GDH=numeric(), tree.id=character(),species=character())


for (i in seq_along(indies)){
  dataoneplant<-filter(HF.bb,tree.id==indies[i])
  burst<-dataoneplant$bb.jd
for(k in c(1:length(burst))){ 
  climatt<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=burst[k]))
  climatt$tree.id<-indies[i]
  climatt$species<-unique(dataoneplant$species)}
df<-rbind(climatt,df) }
  




HF.flo<-filter(HF,!is.na(fopn.jd))
indies2<-unique(HF.flo$tree.id)
df2<-data.frame(Season=character(),End_year=numeric(),Season_days=numeric(),Data_days=numeric(), Perc_complete=numeric(),Chilling_Hours=numeric(),Utah_Model=numeric(),Chill_portions=numeric(), GDH=numeric(), tree.id=character(),species=character())

for (i in seq_along(indies2)){
  dataoneplant<-filter(HF.flo,tree.id==indies2[i])
  burst<-dataoneplant$fopn.jd
  for(k in c(1:length(burst))){ 
    climatt<-as.data.frame(chilling(hourtemps,Start_JDay=1,End_JDay=burst[k]))
    climatt$tree.id<-indies2[i]
    climatt$species<-unique(dataoneplant$species)}
  df2<-rbind(climatt,df2) }

df<-select(df,Season,End_year,GDH,tree.id,species)
df2<-select(df2,Season,End_year,GDH,tree.id,species)

colnames(df)<-c("Season","End_year","GDH_bb","tree.id","species")
colnames(df2)<-c("Season","End_year","GDH_fl","tree.id","species")
daters<-left_join(df2,df)


daters$GDH_diff<-daters$GDH_bb-daters$GDH_fl
?separate()
daters<-separate(daters,tree.id,c("Sp","tree.num"),sep="-")
write.csv(daters,"GDH_diffs_HF.csv",row.names = FALSE)

daters<-read.csv(file = "GDH_diffs_HF.csv",header=TRUE)

sps<-c("FRAM", "ACSA","POTR", "QUVE" ,"QURU" ,"BEPO", "BEPA", "BELE" ,"BEAL" ,"AMSP", "ACRU")

daters$GDD_diff<-daters$GDH_diff/24
daters.hyst<-filter(daters,species %in% sps)
ggplot(daters.hyst,aes(End_year,GDD_diff))+geom_point()+facet_grid(species~tree.num,scale="free")

ggplot(daters.ser,aes(End_year,GDH_diff))+stat_summary(color="blue")+facet_wrap(~species,scale="free")


ggplot(daters.ser,aes(End_year,GDH_diff))+geom_bar(stat="identity",aes(fill=species),position="dodge")+facet_wrap(~species)

ggplot(daters,aes(End_year,GDH_diff))+geom_bar(stat="identity",aes(fill=species),position="dodge")+facet_grid(tree.num~species)

ggplot(daters,aes(End_year,GDH_diff))+geom_bar(stat="identity",aes(fill=species))+facet_wrap(~species)

fullrecords<-c("ACPE","ACRU","ACSA","QURU","BEAL","FAGR","FRAM")
daters.fullrecords<-dplyr::filter(daters, species %in% c(fullrecords))

get_prior(GDD_diff ~ as.factor(End_year)+(as.factor(End_year)|species),data=daters.fullrecords)

daters.fullrecords$End_year<-as.factor(daters.fullrecords$End_year)
mod.fullrec.gdd<-brm(GDH_diff ~ End_year+(End_year|species),
             data=daters.fullrecords,iter=9000,warmup=8000, control = list(adapt_delta=0.99))

modelhere<-mod.fullrec.gdd


###try jsut one species
betal<-filter(daters.fullrecords, species %in% c("BEAL"))
mod.Betal<-brm(GDH_diff ~ End_year+(End_year|tree.num),
                     data=betal,iter=9000,warmup=8000, control = list(adapt_delta=0.99))

Acerub<-filter(daters.fullrecords, species %in% c("ACRU"))
mod.Acru<-brm(GDH_diff ~ End_year,
                     data=Acerub,iter=9000,warmup=8000, control = list(adapt_delta=0.99))

Acepen<-filter(daters.fullrecords, species %in% c("ACPE"))
mod.Acpe<-brm(GDH_diff ~ End_year+End_year:as.factor(tree.num),
              data=Acepen,iter=9000,warmup=8000, control = list(adapt_delta=0.99))




extract_coefs<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.25,0.75,0.975))),"year")
}



beal<-extract_coefs(mod.Betal)
acrb<-extract_coefs(mod.Acru)

new.data.AR<-data.frame(End_year=rep(unique(Acerub$End_year)))#,tree.num=rep(unique(Acerub$tree.num),14))
pred.day.AR<-predict(mod.Acru,newdata=new.data.AR,probs = c(0.1,0.9))
checkpred.AR<-cbind(pred.day.AR,new.data.AR)

pd<-position_dodge(width = 0.8)
ggplot(checkpred.AR,aes(End_year,Estimate))+geom_point(position=pd,color="red")+geom_errorbar(aes(min=Q10,max=Q90),width=0,color="red",position=pd)


tree1<-data.frame(year=1989:2002,FLSdiff=rnorm(14,20,8),GDD.diff=rnorm(14,130,5),tree.id="one")
tree2<-data.frame(year=1989:2002,FLSdiff=rnorm(14,15,10),GDD.diff=rnorm(14,150,5),tree.id="two")
tree3<-data.frame(year=1989:2002,FLSdiff=rnorm(14,15,4),GDD.diff=rnorm(14,160,5),tree.id="three")

fake.trees<-rbind(tree1,tree2,tree3)


scaleFactor <-max(fake.trees$FLSdiff)/ .max(fake.trees$GDD.diff)

a<-ggplot()+stat_summary(aes(fake.trees$year,fake.trees$FLSdiff),color="darkgrey")+
  stat_summary(aes(fake.trees$year,fake.trees$GDD.diff*scaleFactor),color="black")+
  scale_y_continuous("Days between", sec.axis=sec_axis(~./scaleFactor, name="GDD between"))+
  scale_x_discrete("year")+
  theme_bw()+
  theme(
    axis.title.y.left=element_text(color="darkgrey"),
    axis.text.y.left=element_text(color="darkgrey"),
    axis.title.y.right=element_text(color="black"),
    axis.text.y.right=element_text(color="black"))
  

tree1<-data.frame(year=1989:2002,FLSdiff=rnorm(14,20,8),GDD.diff=rnorm(14,130,89),tree.id="one")
tree2<-data.frame(year=1989:2002,FLSdiff=rnorm(14,15,10),GDD.diff=rnorm(14,150,80),tree.id="two")
tree3<-data.frame(year=1989:2002,FLSdiff=rnorm(14,15,4),GDD.diff=rnorm(14,160,100),tree.id="three")

fake.trees2<-rbind(tree1,tree2,tree3)

scaleFactor <-max(fake.trees2$FLSdiff)/max(fake.trees2$GDD.diff)
b<-ggplot()+stat_summary(aes(fake.trees2$year,fake.trees2$FLSdiff),color="darkgrey")+
  stat_summary(aes(fake.trees2$year,fake.trees2$GDD.diff*scaleFactor),color="black")+
  scale_y_continuous("Days between", sec.axis=sec_axis(~./scaleFactor, name="GDD between"))+
  scale_x_discrete("year")+
  theme_bw()+
  theme(
    axis.title.y.left=element_text(color="darkgrey"),
    axis.text.y.left=element_text(color="darkgrey"),
    axis.title.y.right=element_text(color="black"),
    axis.text.y.right=element_text(color="black"))


concept<-ggpubr::ggarrange(a,b)

###real
acerreal<-filter(HF,species=="ACRU")
acerreal<-filter(acerreal,year<=2002)
acerreal$FLS<-acerreal$bb.jd-acerreal$fopn.jd

scaleFactor <-max(acerreal$FLS)/ max(Acerub$GDD_diff)

c<-ggplot()+stat_summary(aes(as.factor(acerreal$year),acerreal$FLS),color="blue")+
  stat_summary(aes(Acerub$End_year,Acerub$GDD_diff*scaleFactor),color="red")+
  scale_y_continuous("Days between", sec.axis=sec_axis(~./scaleFactor, name="GDD between"))+
  scale_x_discrete("year")+
  theme_bw()+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))

jpeg("..//FLOBUDS/Plots/hypothesis1_acerub.jpeg",height=6,width=6,units = "in",res=250)
ggpubr::ggarrange(concept,c,nrow=2)
dev.off()  

new.data.BA<-data.frame(End_year=rep(unique(betal$End_year)),tree.num=rep(unique(betal$tree.num),14))
pred.day.BA<-predict(mod.Betal,newdata=new.data.BA,probs = c(0.1,0.9))
checkpred.BA<-cbind(pred.day.BA,new.data.BA)

ggplot(checkpred.BA,aes(End_year,Estimate))+geom_point(aes(),color="goldenrod1",position)+geom_errorbar(aes(min=Q10,max=Q90),width=0,color="goldenrod1")+facet_wrap(~as.factor(tree.num))

new.data.AP<-data.frame(End_year=rep(unique(Acepen$End_year)),tree.num=rep(unique(Acepen$tree.num),14))
pred.day.AP<-predict(mod.Acpe,newdata=new.data.AP,probs = c(0.1,0.9))
checkpred.AP<-cbind(pred.day.AP,new.data.AP)

ggplot(checkpred.AP,aes(End_year,Estimate))+geom_point(aes(),color="darkgreen")+geom_errorbar(aes(min=Q10,max=Q90),width=0,color="darkgreen")
