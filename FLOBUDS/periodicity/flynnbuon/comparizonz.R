## Dan flynn vs Dan B
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
options(mc.cores = parallel::detectCores())

library(dplyr)
library(brms)

library(dplyr)
library(tidyr)
library(tibble)
library(ggstance)
library(ggthemes)
library(lme4)
library(bayesplot)
library(tidybayes)


setwd("~/Documents/git/proterant/FLOBUDS/periodicity/flynnbuon/")

df<-read.csv("hf314-01-budburst.csv")
db<-read.csv("flobuds_KNB.csv")

###woah df data is needs some cleaning
df<-dplyr::filter(df,site=="HF") #subset to HF data only

### make a common column for GEN.SPA
n <- 4
df$GEN.SPA<-paste(substr(df$sp, 1, n-1), ".", substr(df$sp, n, nchar(df$sp)), sep = "") ## make a column for GEN.SA

ComSp<-intersect(unique(df$GEN.SPA),unique(db$GEN.SPA)) ## select overlapping species

##subset both dataset to species overlap
df<-dplyr::filter(df,GEN.SPA %in% ComSp)
db<-dplyr::filter(db,GEN.SPA %in% ComSp)

##filter flynn to budburst and leadout only (9, 15)
df<-filter(df,tleaf %in% c(6)) ## this is recoded as ``4", see data archeive bud shows green tips

###equalize predictors
df$FORCE<-ifelse(df$warm=="cool",0,10)
df$PHOTO<-ifelse(df$photo=="short",0,4)

db$FORCE<-ifelse(db$Force=="C",0,6)
db$PHOTO<-ifelse(db$Light=="S",0,4)
colnames(db)[5]<-"dayuse" ##leaf exansion bb ch 11


#df<-filter(df,chill=="chill0")
#db<-filter(db,Chill==0)

###first lets see what happens if we just dont account for chilling
df<-dplyr::select(df,GEN.SPA,dayuse,FORCE,PHOTO)
db<-dplyr::select(db,GEN.SPA,dayuse,FORCE,PHOTO)

##dummy variable study
df$study<-0
db$study<-1

dat<-rbind(df,db)
unique(dat$GEN.SPA)


mod<-brm(dayuse~PHOTO+FORCE+study+PHOTO:study+FORCE:study+(1|GEN.SPA),data=dat)
library(xtable)
xtable(fixef(mod,probs = c(.025,.25,.75,.975)))

mod2<-brm(dayuse~PHOTO*FORCE*study+(1|GEN.SPA),data=dat)

fixef(mod2,probs = c(.025,.25,.75,.975))

get_variables(mod)
output<-mod %>%
  spread_draws(b_PHOTO,b_FORCE,`b_PHOTO:study`,`b_FORCE:study`)
###transfor variable to allow interaction ie maineffect = interaction = estimate
output$`b_PHOTO:study`<-output$b_PHOTO+output$`b_PHOTO:study`
output$`b_FORCE:study`<-output$b_FORCE+output$`b_FORCE:study`

output <-output %>% gather("var","estimate",4:7)
output$variable<-ifelse(grepl("PHOTO", output$var),"photoperiod","forcing")  
output$study<-ifelse(grepl("study", output$var),"uncoupled","coupled")

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/modelcomps.jpeg",width = 7,height=8,unit="in",res=300)
ggplot(output,aes(y = variable, x = estimate)) +
  stat_halfeye(aes(fill=study,color=study),alpha=.7)+
  geom_vline(xintercept=0,linetype="dotted")+
  scale_fill_viridis_d()+scale_color_viridis_d()+ggthemes::theme_few()

dev.off() 
#### do it again with the interactions

get_variables(mod2)
output<-mod2 %>%
  spread_draws(b_PHOTO,b_FORCE,`b_PHOTO:FORCE`,`b_PHOTO:study`,`b_FORCE:study`,`b_PHOTO:FORCE:study`)
###transfor variable to allow interaction ie maineffect = interaction = estimate
output$`b_PHOTO:study`<-output$b_PHOTO+output$`b_PHOTO:study`
output$`b_FORCE:study`<-output$b_FORCE+output$`b_FORCE:study`
output$`b_PHOTO:FORCE:study`<-`output$b_PHOTO:FORCE`+output$`b_PHOTO:FORCE:study`

output <-output %>% gather("var","estimate",4:9)
output$variable<-ifelse(grepl("PHOTO", output$var),"photoperiod","forcing")
output$variable<-ifelse(grepl("PHOTO:FORCE", output$var),"interaction",output$variable)

output$study<-ifelse(grepl("study", output$var),"uncoupled","coupled")

jpeg("~/Documents/git/proterant/FLOBUDS/Plots/periodicity_figures/modelcomps_winter.jpeg",width = 7,height=8,unit="in",res=300)
ggplot(output,aes(y = variable, x = estimate)) +
  stat_halfeye(aes(group=study,shape=study),.width=c(.9),alpha=0.75,size=4)+scale_shape_manual(values=c(19,0))+
  geom_vline(xintercept=0,linetype="dotted")+scale_y_discrete(limits = c("interaction", "photoperiod", "forcing"))+
  scale_fill_viridis_d(begin = .4,end = .9)+scale_color_viridis_d(begin = .35,end = .85)+ggthemes::theme_few()+xlab("Phenological sensitivity")+ylab("")+
  xlim(-6,3)

dev.off()





