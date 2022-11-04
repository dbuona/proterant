### Started by Dan Feb 9 2021
### Part 1:### modeling FLS for FNA with measurement models in brms
#FLS ~ pdsi + flower traits + fruit traits
###Part 2: Focusing on prunocerasus with stan model
#based on 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(brms)



#-------------------------#
#------Part 1--------------#
#--------------------------#
setwd("~/Documents/git/proterant/investment/Input/input_clean/")

fruity<-read.csv("fruitsize_clean.csv")
sps<-unique(fruity$specificEpithet)
flowy<-read.csv("petal_clean.csv")
pdsi<-read.csv("pruno_clean_pdsi_wint.csv")
pdsi<-dplyr::filter(pdsi,specificEpithet %in% sps)





fruit.mod<-brm(fruit_diam_mm~1+(1|specificEpithet),data=fruity)
coef(fruit.mod)

fruit.mod<-brm(fruit_diam_mm~1+(1|specificEpithet),data=fruity)
flo.mod<-brm(pental_lengh_mm~1+(1|specificEpithet),data=flowy)
pdsi.mod<-brm(pdsi~1+(1|specificEpithet),data=pdsi)
cold.mod<-brm(wintert~1+(1|specificEpithet),data=pdsi)

flo.out<-as.data.frame(coef(flo.mod))[c(1,2)]
fruit.out<-as.data.frame(coef(fruit.mod))[c(1,2)]
pdsi.out<-as.data.frame(coef(pdsi.mod))[c(1,2)]
cold.out<-as.data.frame(coef(cold.mod))[c(1,2)]

colnames(flo.out)<-c("meanflo","seflo")
colnames(fruit.out)<-c("meanfruit","sefruit")
colnames(pdsi.out)<-c("meanpdsi","sepdsi")
colnames(cold.out)<-c("meancold","secold")

predictor<-cbind(flo.out,fruit.out,pdsi.out,cold.out)
predictor$specificEpithet<-rownames(predictor)

#FLSdata
FLS<-read.csv("FLS_clean.csv")
FLS<-left_join(FLS,predictor)

zscore <- function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}

FLS$flo.z<-zscore(FLS$meanflo)
FLS$fruit.z<-zscore(FLS$meanfruit)
FLS$pdsi.z<-zscore(FLS$meanpdsi)
FLS$cold.z<-zscore(FLS$meancold)
FLS$doy.z<-zscore(FLS$doy)
prunomean.ordz<-brm(bbch.v.scale~doy.z+flo.z+fruit.z+pdsi.z+(1|specificEpithet),
                   data=FLS,
                   family=cumulative("logit"),warmup=3000,iter=4000)

fixef(prunomean.ordz, probs = c(.25,.75))



prunomean.3cats<-brm(bbch.short~doy+flo.z+fruit.z+cold.z+pdsi.z+(1|specificEpithet),
                    data=FLS,
                    family=cumulative("logit"),warmup=3000,iter=4000)

prunomean.gaus<-brm(bbch.v.scale~doy+flo.z+fruit.z+cold.z+pdsi.z+(1|specificEpithet),
                     data=FLS,warmup=3000,iter=4000)
coef(prunomean.gaus,probs = c(.25,.75))
fixef(prunomean.gaus,probs = c(.25,.75))
fixef(prunomean.ordz,probs = c(.25,.75))
fixef(prunomean.3cats,probs = c(.25,.75))

p11<-plot(conditional_effects(prunomean.ordz, "pdsi.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p55<-plot(conditional_effects(prunomean.ordz, "cold.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p22<-plot(conditional_effects(prunomean.ordz, "fruit.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p44<-plot(conditional_effects(prunomean.ordz, "flo.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))

p1<-plot(conditional_effects(prunomean.3cats, "pdsi.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p5<-plot(conditional_effects(prunomean.3cats, "cold.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p2<-plot(conditional_effects(prunomean.3cats, "fruit.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))
p4<-plot(conditional_effects(prunomean.3cats, "flo.z", categorical = TRUE,probs = c(.25,.75),plot=FALSE))


p11<-p11[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p22<-p22[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p44<-p44[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))
p55<-p55[[1]]+ggthemes::theme_few()+scale_color_manual(values=c("hotpink","orange","lightgreen","darkgreen"))+scale_fill_manual(values=c("hotpink","orange","lightgreen","darkgreen"))