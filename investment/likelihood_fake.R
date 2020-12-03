rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(ggplot2)

set.seed(500)
flowers<-c(55,60,63,65,67)
leaves<-c(0,9,11,15,17,19)

species<-(1:15)
#https://stackoverflow.com/questions/24845909/generate-n-random-integers-that-sum-to-m-in-r




d<-data.frame(species=numeric(),bbch.v=numeric(),bbch.f=numeric())

for (i in species){
bbch.f<-sample(flowers,500,replace=TRUE,prob=c(.1,.1,.2,.3,.3))
leafprob<-replicate(1, diff(c(0, sort(runif(5)), 1)))
print(leafprob)
bbch.v<-sample(leaves,500,replace=TRUE,prob=leafprob)
dhere<-data.frame(species=species[i],bbch.v=bbch.v,bbch.f=bbch.f)
d<-rbind(d,dhere)
}

d$FLS<-ifelse(d$bbch.f<67 & d$bbch.v<17,1,0)
d$FLS2<-ifelse(d$bbch.f<67 & d$bbch.v<11,1,0)
###option 1
library(brms)
mod1<-brm(FLS~(1|species),  data=d, family = bernoulli(link = "logit"), warmup = 2000,  iter = 3000, 
          chains = 4)
summary(mod1)
new.data<-data.frame(species=1:15)

mod3<-brm(FLS2~(1|species),  data=d, family = bernoulli(link = "logit"), warmup = 2000,  iter = 3000, 
          chains = 4)

predy3<-fitted(mod3,newdata=new.data)
predy3<-cbind(new.data,predy3)
c<-ggplot(predy3,aes(as.factor(species),Estimate))+geom_bar(stat = "identity")+geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5))+theme_bw()

predy<-fitted(mod1,newdata=new.data)
predy<-cbind(new.data,predy)
a<-ggplot(predy,aes(as.factor(species),Estimate))+geom_bar(stat = "identity")+geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5))+theme_bw()

  mod2<-brm(bbch.v~bbch.f+(bbch.f|species),  data=d, warmup = 2000,  iter = 3000, control=list(adapt_delta=0.95),
          chains = 4)                                                                                         

new.data2<-data.frame(species=1:15,bbch.f=rep(65,15))

predy2<-fitted(mod2,newdata=new.data2)
predy2<-cbind(new.data2,predy2)
b<-ggplot(predy2,aes(as.factor(species),Estimate))+geom_bar(stat = "identity")+geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5))+theme_bw()

ggpubr::ggarrange(a,c,b)
 save.image("likelihoodfake.Rda")




