rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(brms)
library(dplyr)
library(broom)
library(ggplot2)
library(tibble)
library(ggstance)


##1 try changin the exp temperatures
set.seed(1000)
setwd("~/Documents/git/proterant/FLOBUDS")
#load("simulations.Rda")
weather<-data.frame(airt=rep(c(30,10),each=100),doe=rep(1:100,2)) ## set the temperature for each day of the experiment
weather$heatsum<-weather$airt-5 ##this calculate the GDD of each day

weather<-weather %>%group_by(airt) %>% mutate(GDD = cumsum(heatsum)) ###totals the numer of GDD for the whoel experimetn under each treatment 
inds<-1:600 ### now we put 600 cuttings in tyhe chambers

df<-data.frame(tree.id=numeric(),chill=numeric(),photo=numeric(),flowering=numeric(),leafing=numeric()) ## make a data frame
##below simulates the f* for flowering and leaves respective
for(j in c(1:length(inds))){
  photo<-sample(c(0,1),1) ##randomly select low or high photo
  chill<-sample(c(0,1),1)#randomly select low or high chill
  dof<-rnorm(1,200,25)-(rnorm(1,100,20)*chill)-(rnorm(1,20,10)*photo) #notice all that is different from below is F*
  dol<-rnorm(1,400,25)-(rnorm(1,100,20)*chill)-(rnorm(1,20,10)*photo)
  dfhere<-data.frame(tree.id=inds[j],chill=chill,photo=photo,GDD.flo=dof,GDD.leaf=dol)
  
  df<-rbind(df,dfhere)} ##now toyu have the fstars for each indivudal
df$airt<-sample(c(30,10),nrow(df),replace=TRUE)### now assign each a temperature treatment 
weather1<-left_join(weather,df) ## put all you data together

weather.flow<-weather1%>% group_by(tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) ## this pulls out the first time the GDD excede f&
weather.leaf<-weather1%>% group_by(tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) ## this pulls out the first time the GDD excede f&

weather.flow$force<-ifelse(weather.flow$airt==30,1,0)
weather.leaf$force<-ifelse(weather.leaf$airt==30,1,0)

flo.mod.ph<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.flow) # second phase (leaf should ahve higher Temp sense) chilling should be same)
leaf.mod.ph<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.leaf)

extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}

flonodiff<-extract_coefs4HF(flo.mod.ph)
flonodiff$phase<-"flower"
leafnodif<-extract_coefs4HF(leaf.mod.ph)
leafnodif$phase<-"leaf"
nodiff<-rbind(flonodiff,leafnodif)
nodiff<-filter(nodiff,trait!="Intercept")
pd=position_dodgev(height=0.1)
nodiff %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","chill:photo","force:chill","chill","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_rect(ymin=5.5,ymax=6.5,xmin=-80,xmax=-20,fill="hotpink")+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
  ggthemes::theme_few(base_size = 11)+geom_vline(aes(xintercept=0),color="black",linetype="dashed")+
  ylab("cue")+xlab("sensitivity")+ggtitle("Forcing hierarchy" )

##still double


##2 try cahnging the F*s 

df<-data.frame(tree.id=numeric(),chill=numeric(),photo=numeric(),flowering=numeric(),leafing=numeric()) ## make a data frame
##below simulates the f* for flowering and leaves respective
for(j in c(1:length(inds))){
  photo<-sample(c(0,1),1) ##randomly select low or high photo
  chill<-sample(c(0,1),1)#randomly select low or high chill
  dof<-rnorm(1,200,25)-(rnorm(1,100,20)*chill)-(rnorm(1,20,10)*photo) #notice all that is different from below is F*
  dol<-rnorm(1,500,25)-(rnorm(1,100,20)*chill)-(rnorm(1,20,10)*photo)
  dfhere<-data.frame(tree.id=inds[j],chill=chill,photo=photo,GDD.flo=dof,GDD.leaf=dol)
  
  df<-rbind(df,dfhere)} ##now toyu have the fstars for each indivudal
df$airt<-sample(c(30,10),nrow(df),replace=TRUE)### now assign each a temperature treatment 
weather1<-left_join(weather,df) ## put all you data together

weather.flow<-weather1%>% group_by(tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) ## this pulls out the first time the GDD excede f&
weather.leaf<-weather1%>% group_by(tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) ## this pulls out the first time the GDD excede f&

weather.flow$force<-ifelse(weather.flow$airt==30,1,0)
weather.leaf$force<-ifelse(weather.leaf$airt==30,1,0)

flo.mod.ph2<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.flow) # second phase (leaf should ahve higher Temp sense) chilling should be same)
leaf.mod.ph2<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.leaf)


flonodiff2<-extract_coefs4HF(flo.mod.ph2)
flonodiff2$phase<-"flower"
leafnodif2<-extract_coefs4HF(leaf.mod.ph2)
leafnodif2$phase<-"leaf"
nodiff2<-rbind(flonodiff2,leafnodif2)
nodiff2<-filter(nodiff2,trait!="Intercept")
pd=position_dodgev(height=0.1)
nodiff2 %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","chill:photo","force:chill","chill","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_rect(ymin=5.5,ymax=6.5,xmin=-80,xmax=-20,fill="hotpink")+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
  ggthemes::theme_few(base_size = 11)+geom_vline(aes(xintercept=0),color="black",linetype="dashed")+
  ylab("cue")+xlab("sensitivity")+ggtitle("Forcing hierarchy" )
