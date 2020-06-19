rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(brms)
library(dplyr)
library(broom)

setwd("~/Documents/git/proterant/FLOBUDS")
Temp<-25
Tb<-5
threshes<-c(300)
sigma=0.01

days<-1:100

rep<-1:100
chill<-1:2

df<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
    threshhold<-ifelse(chill[j]==1,300+sample(c(20,0,-20),1),300+sample(c(20,0,-20),1)-100)
    for (i in c(1:length(days))){
      y<-ifelse((Temp-Tb)*days[i]==threshhold,days,NA)
      dfhere<-data.frame(temp=25,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
      df<-rbind(df,dfhere)
      df<-df[complete.cases(df),]
      
    }}}
df<-dplyr::select(df,-y)

Temp2<-15 

df2<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
    threshhold<-ifelse(chill[j]==1,300+sample(c(20,0,-20),1),300+sample(c(20,0,-20),1)-100)
    for (i in c(1:length(days))){
      y<-ifelse((Temp2-Tb)*days[i]==threshhold,days,NA)
      df2here<-data.frame(temp=Temp2,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
      df2<-rbind(df2,df2here)
      df2<-df2[complete.cases(df2),]
      
    }}}
df2<-dplyr::select(df2,-y)

df.flo<-rbind(df,df2)
df.flo$phase<-"flowering"


df<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
    threshhold<-ifelse(chill[j]==1,600+sample(c(20,0,-20),1),600+sample(c(20,0,-20),1)-100)
    for (i in c(1:length(days))){
      y<-ifelse((Temp-Tb)*days[i]==threshhold,days,NA)
      dfhere<-data.frame(temp=25,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
      df<-rbind(df,dfhere)
      df<-df[complete.cases(df),]
      
    }}}
df<-dplyr::select(df,-y)

Temp2<-15 

df2<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
    threshhold<-ifelse(chill[j]==1,600+sample(c(20,0,-20),1),600+sample(c(20,0,-20),1)-100)
    for (i in c(1:length(days))){
      y<-ifelse((Temp2-Tb)*days[i]==threshhold,days,NA)
      df2here<-data.frame(temp=Temp2,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
      df2<-rbind(df2,df2here)
      df2<-df2[complete.cases(df2),]
      
    }}}
df2<-dplyr::select(df2,-y)

df.leaf<-rbind(df,df2)
df.leaf$phase<-"leafing"


full<-rbind(df.flo,df.leaf)



library(ggplot2)

scaleFactor <-max(full$days)/ max(full$GDD)

ggplot()+geom_point(data=full,aes(as.factor(temp),days,shape=phase),color="blue")+
  geom_smooth(data=full,method="lm",aes(as.factor(temp),days,group=phase))+
  geom_point(data=full,aes(as.factor(temp),GDD*scaleFactor,shape=phase),color="red")+
  geom_smooth(data=full,method="lm",aes(as.factor(temp),GDD*scaleFactor,group=phase),color="red")+
  scale_y_continuous("Days", sec.axis=sec_axis(~./scaleFactor, name="GDD"))+
  scale_x_discrete("forcing temperature")+
  theme_bw()+facet_wrap(~chill)+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))+ggtitle("Same response")



#####addd chilling
Temp<-25
Tb<-5
threshes<-c(300)
sigma=0.01

days<-1:100

rep<-1:100
chill<-1:2

df<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
  threshhold<-ifelse(chill[j]==1,300+sample(c(20,0,-20),1),300+sample(c(20,0,-20),1)-200)
 for (i in c(1:length(days))){
    y<-ifelse((Temp-Tb)*days[i]==threshhold,days,NA)
    dfhere<-data.frame(temp=25,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
    df<-rbind(df,dfhere)
    df<-df[complete.cases(df),]
    
  }}}
df<-dplyr::select(df,-y)

Temp2<-15 

df2<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

  for (j in chill){
    for (k in c(1:length(rep))){ 
      threshhold<-ifelse(chill[j]==1,300+sample(c(20,0,-20),1),300+sample(c(20,0,-20),1)-300)
  for (i in c(1:length(days))){
    y<-ifelse((Temp2-Tb)*days[i]==threshhold,days,NA)
    df2here<-data.frame(temp=Temp2,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
    df2<-rbind(df2,df2here)
    df2<-df2[complete.cases(df2),]
    
  }}}
df2<-dplyr::select(df2,-y)

df.flo.dif<-rbind(df,df2)
df.flo.dif$phase<-"flowering"


df<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
    threshhold<-ifelse(chill[j]==1,600+sample(c(20,0,-20),1),600+sample(c(20,0,-20),1)-100)
    for (i in c(1:length(days))){
      y<-ifelse((Temp-Tb)*days[i]==threshhold,days,NA)
      dfhere<-data.frame(temp=25,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
      df<-rbind(df,dfhere)
      df<-df[complete.cases(df),]
      
    }}}
df<-dplyr::select(df,-y)

Temp2<-15 

df2<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())

for (j in chill){
  for (k in c(1:length(rep))){ 
    threshhold<-ifelse(chill[j]==1,600+sample(c(20,0,-20),1),600+sample(c(20,0,-20),1)-100)
    for (i in c(1:length(days))){
      y<-ifelse((Temp2-Tb)*days[i]==threshhold,days,NA)
      df2here<-data.frame(temp=Temp2,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
      df2<-rbind(df2,df2here)
      df2<-df2[complete.cases(df2),]
      
    }}}
df2<-dplyr::select(df2,-y)

df.leaf.dif<-rbind(df,df2)
df.leaf.dif$phase<-"leafing"

full1<-rbind(df.leaf.dif,df.flo.dif)



#################  
#full is same resposne to chilling
### full1 is differential respnse


flo.mod.nodif<-brm(days~temp+chill,data=df.flo) # second phase (leaf should ahve higher Temp sense) chilling should be same)
leaf.mod.nodif<-brm(days~temp+chill,data=df.leaf)

flo.mod.dif<-brm(days~temp+chill,data=df.flo.dif)
leaf.mod.dif<-brm(days~temp+chill,data=df.leaf.dif)

extract_coefs4HF<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.025,0.1,0.9,0.975))),"trait")
}

flodiff<-extract_coefs4HF(flo.mod.dif)
flodiff$phase<-"first"
leafdif<-extract_coefs4HF(leaf.mod.dif)
leafdif$phase<-"second"
diff<-rbind(flodiff,leafdif)


flonodiff<-extract_coefs4HF(flo.mod.nodif)
flonodiff$phase<-"first"
leafnodif<-extract_coefs4HF(leaf.mod.nodif)
leafnodif$phase<-"second"
nodiff<-rbind(flonodiff,leafnodif)
diff$model<-"DS"
nodiff$model<-"PH"

simplots<-rbind(nodiff,diff)


simplots<-filter(simplots,trait!="Intercept")


####plot with ggplot
pd=position_dodgev(height=0.4)


  ggplot(simplots,aes(Estimate,trait))+geom_point(aes(Estimate,trait,color=phase),position=pd,size=2)+
  geom_errorbarh(aes(xmin=Q2.5,xmax=Q97.5,color=phase),position=pd,height=0,linetype="dotted")+
  geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),position=pd,height=0,linetype="solid")+
  ggthemes::theme_base(base_size = 10)+geom_vline(aes(xintercept=0),color="black")+facet_wrap(~model)+scale_color_brewer(type = "qual",palette = 7)

#flo.mod.nodif.gdd<-brm(GDD~temp+chill,data=df.flo) #effect on gdd conserved for both phases and cues
#leaf.mod.nodif.gdd<-brm(GDD~temp+chill,data=df.leaf)  

#flo.mod.dif.gdd<-brm(GDD~temp+chill,data=df.flo.dif)
#leaf.mod.dif.gdd-brm(GDD~temp+chill,data=df.leaf.dif)


####plot with plotr:

figpath <- "Plots"

cols <- adjustcolor("indianred3", alpha.f = 0.3) 
my.pal <- rep(brewer.pal(n = 10, name = "Paired"), 8)
# display.brewer.all()
alphahere = 0.4

xlab <- "Model estimate of change in phenophase day"


modelhere <-flo.mod.nodif
modelhere2<-leaf.mod.nodif
source("exp_simcode_brms.R")

modoutput1 <- tidy(modelhere, prob=c(0.5))
modoutput2 <- tidy(modelhere2, prob=c(0.5))

#dev.new()
muplotfx(modelhere,modelhere2, "Precocity hierarchy hypothesis", 8, 4, c(.5,2.5), c(-10, 0))
leg.txt <- c("first phase","second phase")
par(xpd=TRUE) # so I can plot legend outside
legend(1,1,legend=leg.txt,pch=c(19,17),cex=1, bty="n", text.font=3)

dev.off()

###now interceptzoom
muplotfx(modelhere,modelhere2, "intercept zoom1", 8, 8, c(3,4), c(40, 110))

dev.off()
modelhere <-flo.mod.dif
modelhere2<-leaf.mod.dif

modoutput1 <- tidy(modelhere, prob=c(0.5))
modoutput2 <- tidy(modelhere2, prob=c(0.5))

muplotfx(modelhere,modelhere2, "differential sensitivity hypothesis", 8, 4, c(0.5,2.5), c(-20, 0))
leg.txt <- c("first phase","second phase")
par(xpd=TRUE) # so I can plot legend outside
legend(2,1,legend=leg.txt,pch=c(19,17),cex=1, bty="n", text.font=3)

dev.off()


muplotfx(modelhere,modelhere2, "intercept zoom2", 8, 8, c(3,4), c(40, 110))

dev.off()






scaleFactor <-max(full1$days)/ max(full1$GDD)

ggplot()+geom_point(data=full1,aes(as.factor(temp),days,shape=phase),color="blue")+
  geom_smooth(data=full1,method="lm",aes(as.factor(temp),days,group=phase))+
  geom_point(data=full1,aes(as.factor(temp),GDD*scaleFactor,shape=phase),color="red")+
  geom_smooth(data=full1,method="lm",aes(as.factor(temp),GDD*scaleFactor,group=phase),color="red")+
  scale_y_continuous("Days", sec.axis=sec_axis(~./scaleFactor, name="GDD"))+
  scale_x_discrete("forcing temperature")+
  theme_bw()+facet_wrap(~chill)+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))+ggtitle("differential sensitivity")










stop("not an error")

####now add differntial effect of chill
  
  Temp<-25
  Tb<-5
  threshes<-c(300)
  sigma=0.01
  
  days<-1:100
  
  rep<-1:100
  chill<-1:2
  
  df<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())
  
  for (j in chill){
    for (k in c(1:length(rep))){ 
      threshhold<-ifelse(chill[j]==1,300+sample(c(20,0,-20),1),300+sample(c(20,0,-20),1)-120)
      for (i in c(1:length(days))){
        y<-ifelse((Temp-Tb)*days[i]==threshhold,days,NA)
        dfhere<-data.frame(temp=25,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
        df<-rbind(df,dfhere)
        df<-df[complete.cases(df),]
        
      }}}
  df<-dplyr::select(df,-y)
  
  Temp2<-15 
  
  df2<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())
  
  for (j in chill){
    for (k in c(1:length(rep))){ 
      threshhold<-ifelse(chill[j]==1,300+sample(c(20,0,-20),1),300+sample(c(20,0,-20),1)-200)
      for (i in c(1:length(days))){
        y<-ifelse((Temp2-Tb)*days[i]==threshhold,days,NA)
        df2here<-data.frame(temp=Temp2,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
        df2<-rbind(df2,df2here)
        df2<-df2[complete.cases(df2),]
        
      }}}
  df2<-dplyr::select(df2,-y)
  
  df.flo<-rbind(df,df2)
  df.flo$phase<-"flowering"
  
  
  df<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())
  
  for (j in chill){
    for (k in c(1:length(rep))){ 
      threshhold<-ifelse(chill[j]==1,600+sample(c(20,0,-20),1),600+sample(c(20,0,-20),1)-60)
      for (i in c(1:length(days))){
        y<-ifelse((Temp-Tb)*days[i]==threshhold,days,NA)
        dfhere<-data.frame(temp=25,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
        df<-rbind(df,dfhere)
        df<-df[complete.cases(df),]
        
      }}}
  df<-dplyr::select(df,-y)
  
  Temp2<-15 
  
  df2<-data.frame(temp=numeric(),chill=numeric(),days=numeric(),GDD=numeric(),y=numeric(),rep=numeric())
  
  for (j in chill){
    for (k in c(1:length(rep))){ 
      threshhold<-ifelse(chill[j]==1,600+sample(c(20,0,-20),1),600+sample(c(20,0,-20),1)+20)
      for (i in c(1:length(days))){
        y<-ifelse((Temp2-Tb)*days[i]==threshhold,days,NA)
        df2here<-data.frame(temp=Temp2,chill=chill[j],days=days[i],GDD=threshhold,y=y,rep=rep[k])
        df2<-rbind(df2,df2here)
        df2<-df2[complete.cases(df2),]
        
      }}}
  df2<-dplyr::select(df2,-y)
  
  df.leaf<-rbind(df,df2)
  df.leaf$phase<-"leafing"
  
  
  full2<-rbind(df.flo,df.leaf)
  
  
  scaleFactor <-max(full2$days)/ max(full2$GDD)
  
  ggplot()+geom_point(data=full2,aes(as.factor(temp),days,shape=phase),color="blue")+
    geom_smooth(data=full2,method="lm",aes(as.factor(temp),days,group=phase))+
    stat_summary(data=full2,aes(as.factor(temp),GDD*scaleFactor,shape=phase),color="red")+
    geom_smooth(data=full2,method="lm",aes(as.factor(temp),GDD*scaleFactor,group=phase),color="red")+
    scale_y_continuous("Days", sec.axis=sec_axis(~./scaleFactor, name="GDD"))+
    scale_x_discrete("forcing temperature")+
    theme_bw()+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))

  
  
  
  
  ####full1<-rbind(df.flo,df.leaf)
  full1.group<- full1 %>% group_by(temp,chill,phase) %>% summarise(mean.days=mean(days))
  full1.group<-tidyr::spread(full1.group,phase,mean.days)
  full1.group$day.diff<-full1.group$leafing-full1.group$flowering
  
  fullx.group<- full1 %>% group_by(temp,chill,phase) %>% summarise(mean.GDD=mean(GDD))
  fullx.group<-tidyr::spread(fullx.group,phase,mean.GDD)
  fullx.group$GDD.diff<-fullx.group$leafing-fullx.group$flowering
  