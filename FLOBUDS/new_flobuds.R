rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(ggplot2)
library(tidyverse)
library("brms")
library(rstan)
library(arm)
library(rstanarm)
library(tibble)
library(ggstance)
library(survival)
library(sur)
library(survminer)
library(ggthemes)
library("Hmisc")


setwd("~/Documents/git/proterant/FLOBUDS")
load("new_flobud.mods.Rda")
dat<-read.csv("flobudsdata.use.csv",header = TRUE)

dat$Light<-ifelse(dat$Light=="S",0,1)
dat$Force<-ifelse(dat$Force=="C",0,1)

dat<-filter(dat, !GEN.SPA %in% c("AME.SPP","BET.SPP"))
dat

plot.raw.dat<-gather(dat,phase,DOY,2:4)

plot.raw.dat<-unite(plot.raw.dat,Treatcode,5:7,remove=FALSE)
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_1_0"]<-"WL28"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_1_1"]<-"WL56"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_1_0"]<-"CL28"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_0_0"]<-"CS28"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_0_1"]<-"CS56"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="0_1_1"]<-"CL56"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_0_1"]<-"WS56"
plot.raw.dat$Treatcode[plot.raw.dat$Treatcode=="1_0_0"]<-"WS28"
plot.raw.dat$phase[plot.raw.dat$phase=="budburst.9."]<-"budburst"
plot.raw.dat$phase[plot.raw.dat$phase=="leaf_day.15."]<-"expansion"
plot.raw.dat$phase[plot.raw.dat$phase=="flo_day"]<-"flower"

plot.raw.dat.lo<-filter(plot.raw.dat,phase!="budburst")
plot.raw.dat.bb<-filter(plot.raw.dat,phase!="expansion")
goodsp<-c("ACE.PEN","COM.PER","COR.COR","ILE.MUC","PRU.PEN","VAC.COR")
plot.raw.dat.lo<-filter(plot.raw.dat.lo,GEN.SPA %in% c(goodsp) )
plot.raw.dat.bb<-filter(plot.raw.dat.bb,GEN.SPA %in% c(goodsp) )
##clean up variables


loraw<-ggplot(plot.raw.dat.lo,aes(Treatcode,DOY))+geom_point(aes(color=phase,shape=phase),size=.8)+stat_summary(aes(color=phase,shape=phase))+scale_color_manual(values= c("purple","darkgreen"))+scale_shape_manual(values= c(17,15))+facet_wrap(~GEN.SPA)+theme(axis.text.x = element_text(angle=45))+theme_linedraw()
bbraw<-ggplot(plot.raw.dat.bb,aes(Treatcode,DOY))+geom_point(aes(color=phase,shape=phase),size=.8)+stat_summary(aes(color=phase,shape=phase))+facet_wrap(~GEN.SPA)+scale_color_manual(values= c("purple","darkgreen"))+scale_shape_manual(values= c(17,0))+theme(axis.text.x = element_text(size = .3))+theme_linedraw()
grid.arrange(bbraw,loraw)

plot.raw.dat$taxa<-NA
plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="ACE.PEN"]<-"A. pensylvanicum"
plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="COM.PER"]<-"C. peregrina"
plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="COR.COR"]<-"C. cornuta"
plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="ILE.MUC"]<-"I. mucronata"
plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="PRU.PEN"]<-"P. pensylvanica"
plot.raw.dat$taxa[plot.raw.dat$GEN.SPA=="VAC.COR"]<-"V corymbosum"

spp.raw<-ggplot(plot.raw.dat,aes(Treatcode,DOY))+geom_point(aes(color=phase,shape=phase),size=.8)+stat_summary(aes(color=phase,shape=phase))+scale_shape_manual(values= c(0,2,6))+facet_wrap(~GEN.SPA)+theme(axis.text.x = element_text(angle=45))

goodspfull<-filter(plot.raw.dat,GEN.SPA %in% c(goodsp)) 
goodspfull$force<-ifelse(goodspfull$Force==1,"High","Low")
goodspfull$photo<-ifelse(goodspfull$Light==1,"Long day","Short day")
goodspfull$chill<-ifelse(goodspfull$Chill==1,"High chill","Low chill")

ggplot(goodspfull,aes(taxa,DOY))+geom_point(aes(color=phase,shape=force),size=0.8)+stat_summary(aes(color=phase,shape=force))+facet_grid(chill~photo)+theme_base()+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Treatment combination")+theme_base(base_size = 7)

pdf("Plots/flo_buds_figures/goodsps_rawplots2.pdf",width=11, height=6)
#ggplot(goodspfull,aes(Treatcode,DOY))+geom_point(aes(color=phase,shape=phase),size=.8)+stat_summary(aes(color=phase,shape=phase))+facet_wrap(~taxa)+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Treatment combination")+theme_base(base_size = 7)
#ggplot(goodspfull,aes(force,DOY))+geom_point(aes(color=phase,shape=photo),size=0.8)+stat_summary(aes(color=phase,shape=photo))+facet_grid(chill~taxa)+theme_base()+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Treatment combination")+theme_base(base_size = 12)
ggplot(goodspfull,aes(phase,DOY))+geom_point(aes(color=force,shape=photo),size=0.8)+stat_summary(aes(color=force,shape=photo))+facet_grid(chill~taxa)+theme_base()+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Phenophase")+theme_base(base_size = 8)

#ggplot(goodspfull,aes(force,DOY))+geom_point(aes(color=taxa,shape=phase),size=0.8)+stat_summary(aes(color=taxa,shape=phase))+facet_grid(chill~photo)+theme_base()+theme(axis.text.x = element_text(angle=-45,size=8))+ylab("Day of experiment")+xlab("Treatment combination")+theme_base(base_size = 7)
dev.off()

bb.int<-get_prior(budburst.9.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())
mod.bb.int<-brm(budburst.9. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                   data = dat, family = gaussian(),
                   iter= 4000,
                   warmup = 3000)   
summary(mod.bb.int)

flo.int<-get_prior(flo_day~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())
mod.flo.int<-brm(flo_day~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3200)   
summary(mod.flo.int)

extract_coefs<-function(x){rownames_to_column(as.data.frame(fixef(x, summary=TRUE,probs=c(0.10,.25,.75,0.90))),"Predictor")
}
extract_ranef<-function(x){tibble::rownames_to_column(as.data.frame(ranef(x, summary=TRUE,probs=c(0.10,.25,.75,0.90))),"GEN.SPA")
}


flowy<-extract_coefs(mod.flo.int)
leafy<-extract_coefs(mod.bb.int)
leafy$phase<-"foliate"
flowy$phase<-"floral"

bothy<-rbind(flowy,leafy)
bothy<-filter(bothy,Predictor!="Intercept")
bothy$Predictor[bothy$Predictor=="Force"]<-"main:Force"
bothy$Predictor[bothy$Predictor=="Light"]<-"main:Light"
bothy$Predictor[bothy$Predictor=="Chill"]<-"main:Chill"
bothy$Predictor[bothy$Predictor=="Chill:Light"]<-"int:Chill:Light"
bothy$Predictor[bothy$Predictor=="Light:Force"]<-"int:Light:Force"
bothy$Predictor[bothy$Predictor=="Chill:Force"]<-"int:Chill:Force"


bothy$phase[bothy$phase=="floral"]<-"flower open"
bothy$phase[bothy$phase=="foliate"]<-"leaf budburst"
bothy$phase[bothy$phase=="foliate-lo"]<-"leaf expansion"


pdf("Plots/flo_buds_figures/fobb_maineffects.pdf",width=8.5, height=6)
pd2=position_dodgev(height=0.3)
ggplot(bothy,aes(Estimate,Predictor))+geom_point(aes(color=phase,shape=phase),position=pd2, size=4)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()
##random effects

ran.flo<-extract_ranef(mod.flo.int)
colnames(ran.flo)

a<-dplyr::select(ran.flo,1:7)
b<-dplyr::select(ran.flo,1,8:13)
c<-dplyr::select(ran.flo,1,14:19)
d<-dplyr::select(ran.flo,1,20:25)
e<-dplyr::select(ran.flo,1,26:31)
f<-dplyr::select(ran.flo,1,32:37)
g<-dplyr::select(ran.flo,1,38:43)

a$Predictor<-"Intercept"
b$Predictor<-"main:chilling"
c$Predictor<-"main:light"
d$Predictor<-"main:forcing"
e$Predictor<-"int:chillxlight"
f$Predictor<-"int:chillxforce"
g$Predictor<-"int:lightxforce"

call<-c("GEN.SPA","Estimate","Error","Q10","Q25","Q75","Q90","Predictor")
colnames(a)<-call
colnames(b)<-call
colnames(c)<-call
colnames(d)<-call
colnames(e)<-call
colnames(f)<-call
colnames(g)<-call

a[,2:7]<-a[,2:7]+flowy[1,2]
b[,2:7]<-b[,2:7]+flowy[2,2]
c[,2:7]<-c[,2:7]+flowy[3,2]
d[,2:7]<-d[,2:7]+flowy[4,2]
e[,2:7]<-e[,2:7]+flowy[5,2]
f[,2:7]<-f[,2:7]+flowy[6,2]
g[,2:7]<-g[,2:7]+flowy[7,2]



flo.sps<-rbind(a,b,c,d,e,f,g)
flo.sps$phase<-"floral"
pd2=position_dodgev(height=0.6)
#flo.sps<-filter(flo.sps,Predictor!="xIntercept")

ggplot(flo.sps,aes(Estimate,Predictor))+geom_point(aes(color=GEN.SPA),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=GEN.SPA),linetype="solid",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()

###now leaf estiamte
ran.leaf<-extract_ranef(mod.bb.int)

aa<-dplyr::select(ran.leaf,1:7)
bb<-dplyr::select(ran.leaf,1,8:13)
cc<-dplyr::select(ran.leaf,1,14:19)
dd<-dplyr::select(ran.leaf,1,20:25)
ee<-dplyr::select(ran.leaf,1,26:31)
ff<-dplyr::select(ran.leaf,1,32:37)
gg<-dplyr::select(ran.leaf,1,38:43)

aa$Predictor<-"Intercept"
bb$Predictor<-"main:chilling"
cc$Predictor<-"main:light"
dd$Predictor<-"main:forcing"
ee$Predictor<-"int:chillxlight"
ff$Predictor<-"int:chillxforce"
gg$Predictor<-"int:lightxforce"

call<-c("GEN.SPA","Estimate","Error","Q10","Q25","Q75","Q90","Predictor")
colnames(aa)<-call
colnames(bb)<-call
colnames(cc)<-call
colnames(dd)<-call
colnames(ee)<-call
colnames(ff)<-call
colnames(gg)<-call


aa[,2:7]<-aa[,2:7]+leafy[1,2]
bb[,2:7]<-bb[,2:7]+leafy[2,2]
cc[,2:7]<-cc[,2:7]+leafy[3,2]
dd[,2:7]<-dd[,2:7]+leafy[4,2]
ee[,2:7]<-ee[,2:7]+leafy[5,2]
ff[,2:7]<-ff[,2:7]+leafy[6,2]
gg[,2:7]<-gg[,2:7]+leafy[7,2]

leaf.sps<-rbind(aa,bb,cc,dd,ee,ff,gg)
leaf.sps$phase<-"foliate"
#leaf.sps<-filter(leaf.sps, Predictor!="xIntercept")

ggplot(leaf.sps,aes(Estimate,Predictor))+geom_point(aes(color=GEN.SPA),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=GEN.SPA),linetype="solid",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()
sps.plot<-rbind(leaf.sps,flo.sps)

sps.plot$fls<-ifelse(sps.plot$GEN.SPA %in% c("ACE.RUB","COM.PER","COR.COR"),"hyst","ser")


jpeg("Plots/flo_buds_figures/fobb_sppeffects.wintercept.jpeg")
pd2=position_dodgev(height=.6)
ggplot(sps.plot,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()

#without intecepts
sps.plot.effectonly<-filter(sps.plot, Predictor!="Intercept")
jpeg("Plots/flo_buds_figures/fobb_sppeffects.nointercept.jpeg")
ggplot(sps.plot.effectonly,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()
########Leafout vs flowering
lo.int<-get_prior(leaf_day.15.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian())
mod.lo.int<-brm(leaf_day.15. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3000)   
summary(mod.lo.int)

leafouty<-extract_coefs(mod.lo.int)
leafouty$phase<-"foliate-lo"
flowy$phase<-"floral"

bothy<-rbind(flowy,leafouty,leafy)
bothy<-filter(bothy,Predictor!="Intercept")
bothy$Predictor[bothy$Predictor=="Force"]<-"main:Force"
bothy$Predictor[bothy$Predictor=="Light"]<-"main:Light"
bothy$Predictor[bothy$Predictor=="Chill"]<-"main:Chill"
bothy$Predictor[bothy$Predictor=="Chill:Light"]<-"int:Chill:Light"
bothy$Predictor[bothy$Predictor=="Light:Force"]<-"int:Light:Force"
bothy$Predictor[bothy$Predictor=="Chill:Force"]<-"int:Chill:Force"

jpeg("Plots/flo_buds_figures/folo_maineffects.jpeg")
pd2=position_dodgev(height=0.3)
ggplot(bothy,aes(Estimate,Predictor))+geom_point(aes(color=phase,shape=phase),position=pd2, size=4)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()+scale_color_manual(values=c("darkgrey", "black"))+ggtitle("Main effects flower vs. leafout")
dev.off()

ran.leafout<-extract_ranef(mod.lo.int)

aaa<-dplyr::select(ran.leafout,1:7)
bbb<-dplyr::select(ran.leafout,1,8:13)
ccc<-dplyr::select(ran.leafout,1,14:19)
ddd<-dplyr::select(ran.leafout,1,20:25)
eee<-dplyr::select(ran.leafout,1,26:31)
fff<-dplyr::select(ran.leafout,1,32:37)
ggg<-dplyr::select(ran.leafout,1,38:43)

aaa$Predictor<-"Intercept"
bbb$Predictor<-"main:chilling"
ccc$Predictor<-"main:light"
ddd$Predictor<-"main:forcing"
eee$Predictor<-"int:chillxlight"
fff$Predictor<-"int:chillxforce"
ggg$Predictor<-"int:lightxforce"

call<-c("GEN.SPA","Estimate","Error","Q10","Q25","Q75","Q90","Predictor")
colnames(aaa)<-call
colnames(bbb)<-call
colnames(ccc)<-call
colnames(ddd)<-call
colnames(eee)<-call
colnames(fff)<-call
colnames(ggg)<-call


aaa[,2:7]<-aaa[,2:7]+leafouty[1,2]
bbb[,2:7]<-bbb[,2:7]+leafouty[2,2]
ccc[,2:7]<-ccc[,2:7]+leafouty[3,2]
ddd[,2:7]<-ddd[,2:7]+leafouty[4,2]
eee[,2:7]<-eee[,2:7]+leafouty[5,2]
fff[,2:7]<-fff[,2:7]+leafouty[6,2]
ggg[,2:7]<-ggg[,2:7]+leafouty[7,2]

leafout.sps<-rbind(aaa,bbb,ccc,ddd,eee,fff,ggg)
leafout.sps$phase<-"foliate-lo"

sps.plot2<-rbind(leafout.sps,flo.sps)
jpeg("Plots/flo_buds_figures/folo_sppeffects.wintercept.jpeg")
pd2=position_dodgev(height=.6)
ggplot(sps.plot2,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()


sps.plot3<-rbind(leafout.sps,leaf.sps,flo.sps)
jpeg("Plots/flo_buds_figures/3phase_sppeffects_forsupplimnet.jpeg")
pd2=position_dodgev(height=.6)
ggplot(sps.plot3,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()

goodsps2<-dplyr::filter(sps.plot3, GEN.SPA %in% c("ACE.PEN","COM.PER","COR.COR","ILE.MUC","PRU.PEN","VAC.COR"))

goodsps2$phase[goodsps2$phase=="floral"]<-"flower open"
goodsps2$phase[goodsps2$phase=="foliate"]<-"leaf budburst"
goodsps2$phase[goodsps2$phase=="foliate-lo"]<-"leaf expansion"
goodsps2$taxa<-NA
goodsps2$taxa[goodsps2$GEN.SPA=="ACE.PEN"]<-"A. pensylvanicum"
goodsps2$taxa[goodsps2$GEN.SPA=="COM.PER"]<-"C. peregrina"
goodsps2$taxa[goodsps2$GEN.SPA=="COR.COR"]<-"C. cornuta"
goodsps2$taxa[goodsps2$GEN.SPA=="ILE.MUC"]<-"I. mucronata"
goodsps2$taxa[goodsps2$GEN.SPA=="PRU.PEN"]<-"P. pensylvanica"
goodsps2$taxa[goodsps2$GEN.SPA=="VAC.COR"]<-"V. corymbosum"

jpeg("Plots/flo_buds_figures/3phase_sppeffects_goodsps.jpeg")
pd2=position_dodgev(height=.6)
ggplot(goodsps2,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=3)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~GEN.SPA)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()

sps.plot.3effectonly<-filter(goodsps2, Predictor!="Intercept")

jpeg("Plots/flo_buds_figures/3phase_sppeffectsgoodsps_nointercept.jpeg",height=500)
ggplot(sps.plot.3effectonly,aes(Estimate,Predictor))+geom_point(aes(shape=phase,color=phase),position=pd2, size=2)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.5)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+facet_wrap(~taxa)+geom_vline(aes(xintercept=0),color="black")+theme_base(base_size = 10)
dev.off()

bothy<-rbind(flowy,leafouty,leafy)
bothy<-filter(bothy,Predictor!="Intercept")
bothy$Predictor[bothy$Predictor=="Force"]<-"main:Force"
bothy$Predictor[bothy$Predictor=="Light"]<-"main:Light"
bothy$Predictor[bothy$Predictor=="Chill"]<-"main:Chill"
bothy$Predictor[bothy$Predictor=="Chill:Light"]<-"int:Chill:Light"
bothy$Predictor[bothy$Predictor=="Light:Force"]<-"int:Light:Force"
bothy$Predictor[bothy$Predictor=="Chill:Force"]<-"int:Chill:Force"

jpeg("Plots/flo_buds_figures/2phase_maineffects.jpeg")
pd2=position_dodgev(height=0.3)
ggplot(bothy,aes(Estimate,Predictor))+geom_point(aes(color=phase,shape=phase),position=pd2, size=4)+geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=phase),linetype="solid",position=pd2,width=0,size=0.7)+geom_errorbarh(aes(xmin=Q10,xmax=Q90,color=phase),linetype="dotted",position=pd2,width=0,size=0.7)+geom_vline(aes(xintercept=0),color="black")+theme_base()
dev.off()


####Do hysteranthous species flower earlier than no hysteranthous

FLS<-filter(dat, !GEN.SPA %in% c("ACE.SAC","BET.ALL"))
unique(FLS$GEN.SPA)
hys<-c("ACE.RUB","COM.PER","COR.COR","ACE.SAC","BET.ALL")
syn<-c("ACE.PEN","PRU.PEN","VAC.COR")
ser<-c("ILE.MUC","ILE.VER","PRU.VIR","VIB.ACE")

colnames(dat)
meanf<-dat %>% group_by(GEN.SPA)%>% dplyr::summarise(meano=mean(flo_day,na.rm=TRUE))
meanl<-dat %>% group_by(GEN.SPA)%>% dplyr::summarise(meanl=mean(budburst.9.,na.rm=TRUE))
averages<-left_join(meanf,meanl)
averages$FLS.cont<-averages$meanl-averages$meano
dat<-left_join(dat,averages,by="GEN.SPA")
dat$FLS.bin.phys<-ifelse(dat$GEN.SPA %in% hys,1,0)
colnames(dat)
summary(lm(flo_day~FLS.bin.phys,data=dat))
summary(lm(leaf_day.15.~FLS.bin.phys,data=dat))
summary(lm(budburst.9.~FLS.bin.phys,data=dat))

summary(lm(flo_day~FLS.cont,data=dat))
summary(lm(leaf_day.15.~FLS.cont,data=dat))
summary(lm(budburst.9.~FLS.cont,data=dat))

?summarise()
newdata<-data.frame(Chill = c(0,0,1,1,0,1,1,0),
                      Light = c(0,1,0,1,0,1,0,1),
                      Force=  c(0,1,1,0,1,1,0,0),
                      GEN.SPA=c("VAC.COR"))
                                
newdata<-unite(newdata,treatcode,1:3,remove=FALSE)

bb.pred<-predict(mod.bb.int, newdata = newdata,summary=TRUE,   probs = c(0.025, 0.975))
flo.pred<-predict(mod.flo.int, newdata = newdata,summary=TRUE,   probs = c(0.025, 0.975))



save.image("new_flobud.mods.Rda")
