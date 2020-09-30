rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(brms)
library(dplyr)
library(broom)
library(ggplot2)
library(tibble)
library(ggstance)

set.seed(1000)
setwd("~/Documents/git/proterant/FLOBUDS")
weather<-data.frame(airt=rep(c(25,15),each=100),doe=rep(1:100,2)) ## set the temperature for each day of the experiment
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
df$airt<-sample(c(25,15),nrow(df),replace=TRUE)### now assign each a temperature treatment 
weather1<-left_join(weather,df) ## put all you data together

weather.flow<-weather1%>% group_by(tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) ## this pulls out the first time the GDD excede f&
weather.leaf<-weather1%>% group_by(tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) ## this pulls out the first time the GDD excede f&
##########
df2<-data.frame(tree.id=numeric(),chill=numeric(),photo=numeric(),flowering=numeric(),leafing=numeric())

for(j in c(1:length(inds))){ #now weill make some differntial senstivity
  photo<-sample(c(0,1),1)
  chill<-sample(c(0,1),1)
  dof<-rnorm(1,200,25)-(rnorm(1,100,25)*chill)-(rnorm(1,20,10)*photo) #here fstar are dependent of chillign difference and forcing too
  dol<-rnorm(1,400,25)-(rnorm(1,200,25)*chill)-(rnorm(1,0,10)*photo)
  df2here<-data.frame(tree.id=inds[j],chill=chill,photo=photo,GDD.flo=dof,GDD.leaf=dol)
  
  df2<-rbind(df2,df2here)}
df2$airt<-sample(c(25,15),nrow(df2),replace=TRUE)
weather2<-left_join(weather,df2)

weather.flow.2<-weather2%>% group_by(tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) 
weather.leaf.2<-weather2%>% group_by(tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) 



df3<-data.frame(tree.id=numeric(),chill=numeric(),photo=numeric(),flowering=numeric(),leafing=numeric())

for(j in c(1:length(inds))){ #now weill make some differntial senstivity
  photo<-sample(c(0,1),1)
  chill<-sample(c(0,1),1)
  dof<-rnorm(1,400,25)-(rnorm(1,100,25)*chill)-(rnorm(1,20,10)*photo) #here fstar are dependent of chillign difference and forcing too
  dol<-rnorm(1,400,25)-(rnorm(1,200,25)*chill)-(rnorm(1,0,10)*photo)
  df3here<-data.frame(tree.id=inds[j],chill=chill,photo=photo,GDD.flo=dof,GDD.leaf=dol)
  
  df3<-rbind(df3,df3here)}
df3$airt<-sample(c(25,15),nrow(df3),replace=TRUE)
weather3<-left_join(weather,df3)

weather.flow.3<-weather3%>% group_by(tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) 
weather.leaf.3<-weather3%>% group_by(tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) 


df4<-data.frame(tree.id=numeric(),chill=numeric(),photo=numeric(),flowering=numeric(),leafing=numeric())

for(j in c(1:length(inds))){ #now weill make some differntial senstivity
  photo<-sample(c(0,1),1)
  chill<-sample(c(0,1,2),1)
  dof<-ifelse(chill==2,rnorm(1,200,25)-(rnorm(1,100,25)*chill)-(rnorm(1,20,10)*photo),rnorm(1,200,25)-(rnorm(1,20,10)*photo)) #here fstar are dependent of chillign difference and forcing too
  dol<-ifelse(chill!=0,rnorm(1,400,25)-(rnorm(1,100,25)*chill)-(rnorm(1,0,10)*photo),rnorm(1,400,25)-(rnorm(1,0,10)*photo))
  df4here<-data.frame(tree.id=inds[j],chill=chill,photo=photo,GDD.flo=dof,GDD.leaf=dol)
  
  df4<-rbind(df4,df4here)}
df4$airt<-sample(c(25,15),nrow(df4),replace=TRUE)
weather4<-left_join(weather,df4)

weather.flow.4<-weather4%>% group_by(tree.id) %>% filter(GDD >=GDD.flo) %>% slice(1) 
weather.leaf.4<-weather4%>% group_by(tree.id) %>% filter(GDD >=GDD.leaf) %>% slice(1) 



weather.flow$phase<-"flower"
weather.flow.2$phase<-"flower"
weather.flow.3$phase<-"flower"
weather.flow.4$phase<-"flower"
weather.leaf$phase<-"leaf"
weather.leaf.2$phase<-"leaf"
weather.leaf.3$phase<-"leaf"
weather.leaf.4$phase<-"leaf"


extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}

###you need to dummy code airt before you keep going on this
weather.flow.4$force<-ifelse(weather.flow.4$airt==25,1,0)
weather.flow.3$force<-ifelse(weather.flow.3$airt==25,1,0)
weather.flow.2$force<-ifelse(weather.flow.2$airt==25,1,0)
weather.flow$force<-ifelse(weather.flow$airt==25,1,0)

weather.leaf.4$force<-ifelse(weather.leaf.4$airt==25,1,0)
weather.leaf.3$force<-ifelse(weather.leaf.3$airt==25,1,0)
weather.leaf.2$force<-ifelse(weather.leaf.2$airt==25,1,0)
weather.leaf$force<-ifelse(weather.leaf$airt==25,1,0)

###model the sensitivity


flo.mod.dif.3<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.flow.3)
leaf.mod.dif.3<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.leaf.3)


flo.mod.dif.int<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.flow.2)
leaf.mod.dif.int<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.leaf.2)


flo.mod.ph<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.flow) # second phase (leaf should ahve higher Temp sense) chilling should be same)
leaf.mod.ph<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.leaf)


phh.thresh.flo<-filter(weather.flow.4,chill==2)
phh.thresh.leaf<-filter(weather.leaf.4,chill %in%c(1,2))

flo.mod.thresh<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.flow.4)
leaf.mod.thresh<-brm(doe~force+chill+photo+force:chill+force:photo+photo:chill,data=weather.leaf.4)

flo.mod.thresh2<-brm(doe~force+photo+force:photo,data=phh.thresh.flo)
leaf.mod.thresh2<-brm(doe~force+photo+force:photo,data=phh.thresh.leaf)

flodiff.int<-extract_coefs4HF(flo.mod.dif.int)
flodiff.int$phase<-"flower"
leafdif.int<-extract_coefs4HF(leaf.mod.dif.int)
leafdif.int$phase<-"leaf"
diff.int<-rbind(flodiff.int,leafdif.int)


flothresh<-extract_coefs4HF(flo.mod.thresh)
flothresh$phase<-"flower"
leafthresh<-extract_coefs4HF(leaf.mod.thresh)
leafthresh$phase<-"leaf"
thresh<-rbind(flothresh,leafthresh)

flothresh2<-extract_coefs4HF(flo.mod.thresh2)
flothresh2$phase<-"flower"
leafthresh2<-extract_coefs4HF(leaf.mod.thresh2)
leafthresh2$phase<-"leaf"
thresh2<-rbind(flothresh2,leafthresh2)


flodiff.3<-extract_coefs4HF(flo.mod.dif.3)
flodiff.3$phase<-"flower"
leafdif.3<-extract_coefs4HF(leaf.mod.dif.3)
leafdif.3$phase<-"leaf"
diff.3<-rbind(flodiff.3,leafdif.3)

flonodiff<-extract_coefs4HF(flo.mod.ph)
flonodiff$phase<-"flower"
leafnodif<-extract_coefs4HF(leaf.mod.ph)
leafnodif$phase<-"leaf"
nodiff<-rbind(flonodiff,leafnodif)




nodiff<-filter(nodiff,trait!="Intercept")
diff.int<-filter(diff.int,trait!="Intercept")
diff.3<-filter(diff.3,trait!="Intercept")
thresh<-filter(thresh,trait!="Intercept")
thresh2<-filter(thresh2,trait!="Intercept")

pd=position_dodgev(height=0.1)
a<-nodiff %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","chill:photo","force:chill","chill","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
 ggthemes::theme_clean(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  ggtitle("Precocity hierarchy")+ylab("cue")+xlab("sensitivity")+xlim(-25,15)

b<-diff.int %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","chill:photo","force:chill","chill","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
  ggthemes::theme_clean(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  ggtitle("Precocity hierarchy & differential sensitity")+ylab("cue")+xlab("sensitivity")+xlim(-25,15)


c<-diff.3 %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","chill:photo","force:chill","chill","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
  ggthemes::theme_clean(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  ggtitle("Differential sensitivity")+ylab("cue")+xlab("sensitivity")+xlim(-25,15)

d<-thresh %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","chill:photo","force:chill","chill","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
  ggthemes::theme_clean(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  ggtitle("Differential sensitivity-threshholds" )+ylab("cue")+xlab("sensitivity")

e<-thresh2 %>%
  arrange(Estimate) %>%
  mutate(trait = factor(trait, levels=c("force:photo","photo","force"))) %>%
  ggplot(aes(Estimate,trait))+geom_point(aes(shape=phase),position=pd,size=3)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,group=phase),position=pd,width=0)+
  ggthemes::theme_clean(base_size = 11)+geom_vline(aes(xintercept=0),color="black")+
  ggtitle("Differential sensitivity-threshholds" )+ylab("cue")+xlab("sensitivity")


png("Plots/Flobuds_manuscript_figs/simulations.png",width = 5,height = 5,units = "in",res=200)
ggpubr::ggarrange(a,b,c,ncol=1,common.legend = TRUE,legend="bottom",labels = c("a","b","c","d"))
dev.off()
