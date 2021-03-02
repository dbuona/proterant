## can you just combin dposterions?
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(brms)
library(ggplot2)
   #####################
###joint model thinking ########
##########################
#---------IDEA 1------
#model 1
#FLS ~ 1|species
###model 2
#pdsi ~ 1|species
###joint model
#pdsi~FLS or (pdsi ~ 1|species) ~ (FLS ~ 1|species)

#----IDEA 2------------
#model 1:
#FLS ~1|species

##joint model
#pdsi ~ (1|[FLS~1|species])

#----IDEA 3------
# do we even need joint modeling
# pdsi~FLS+(FLS|species)


set.seed(2021)
##vary intercept


nsp = 14 # number of species
ntot = 100 # numbers of obs per species.

baseinter <- 3 # baseline intercept (budvol) across all species
spint <- baseinter +sample(-2:2,length(1:nsp),replace=T) # different intercepts by species
spintvar
#
  
df<-data.frame(species=numeric(),FLS=numeric(),pdsi=numeric())

for (i in spint){
  pdsi=spint[i]+
    for (i in nsp){
    dfhere<-data.frame(species=sp[i],pdsi=pdsi)
                   df<-rbind(df,dfhere)
  }

ggplot(df,aes(as.factor(species),pdsi))+stat_summary()

df2<-data.frame(species=numeric(),FLS=numeric())
for (i in sp){
    FLS<-sample(1:5,size = 1000,replace = TRUE,prob = sample(c(.1,.25,.05,.1,.7),5,replace=FALSE))
  df2here<-data.frame(species=sp[i],FLS=FLS)
  df2<-rbind(df2,df2here)
}

ggplot(df2,aes(as.factor(species),FLS))+stat_summary()

ggplot()+
  geom_point(data=df,aes(species,pdsi),color="black")+
  geom_point(data=df2,aes(species,FLS-4.5),color="firebrick")+
  scale_y_continuous( sec.axis=sec_axis(~.+4.5, name="FLS"))+
  geom_smooth(data=df,method="lm",aes(species,pdsi),color="black")+
 geom_smooth(data=df2,method="lm",aes(species,FLS-4.5),color="firebrick")

cor(df$pdsi,df2$FLS)
###Okay above reflect data collect, but in truth we should try and generate the relationship between pdsi and FLS more directly
df2$pdsi<--1.5-.18*df2$species+rnorm(nrow(df2),0,1)

ggplot()+
  stat_summary(data=df2,aes(as.factor(species),pdsi),color="black")+
  stat_summary(data=df2,aes(as.factor(species),FLS-4.5),color="firebrick")+
  scale_y_continuous( sec.axis=sec_axis(~.+4.5, name="FLS"))

library(lme4)
mod<-lmer(df2$pdsi~df2$FLS+(1|df2$species))
summary(mod)

#######
stop("below is scratch")
setwd("~/Documents/git/proterant/investment/input")
load("PrunusFLSs.Rda")

goomod<-brm(pdsi~bbch.v+(1|specificEpithet),data=d.flo)
summary(goomod)



sampleFLS <- posterior_samples(mod2a)### extract the posteriors
samplePDSI<-posterior_samples(mod.pdsi) 
## not sure which columns to use
pdsi.goo<-dplyr::select(mod.pdsi.out,species,Estimate, SE)
colnames(pdsi.goo)<-c("species","pdsi","pdsi.se")

FLS.goo<-dplyr::select(mod2a.ord.out,species,Estimate, SE)
colnames(FLS.goo)<-c("species","FLS","FLS.se")
goo.dat<-dplyr::left_join(FLS.goo,pdsi.goo)

goomod<- brm(FLS~ me(pdsi, pdsi.se), data = goo.dat, 
            save_mevars = TRUE)

bform <- bf(me(pdsi, pdsi.se) ~ me(FLS, FLS.se)) + set_mecor(FALSE)
fit2 <- brm(bform, data = goo.dat, save_mevars = TRUE)
summary(fit2)

## can't put this on the right side of equation



