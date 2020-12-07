load("writing/flobud.main.mods.Rda")
small<-filter(dat,!is.na(flo_day))
small<-filter(small,!is.na(budburst.9.))

small$FLS<-small$budburst.9.-small$flo_day

small.fls<-get_prior(FLS~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = small, family = gaussian())
mod.FLS.small<-brm(FLS ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                   data = small, family = gaussian(),control = list(adapt_delta = 0.95),
                   iter= 7000,
                   warmup = 6000)





###prediction plots
HFreal<-read.csv(file = "..//Data/hf003-05-mean-ind.csv")
HFreal$FLS<-HFreal$bb.jd-HFreal$fopn.jd
HFrealmeans<-HFreal %>% group_by(species) %>% dplyr::summarize(Estimate=mean(FLS,na.rm = TRUE),Q94.5=max(FLS,na.rm = TRUE),Q5.5=min(FLS,na.rm = TRUE))
HFrealmeans<-filter(HFrealmeans,species %in% c('ACPE',"ACRU","NEMU","ILVE","VACO"))
HFrealmeans$GEN.SPA<-NA

HFrealmeans$GEN.SPA[which(HFrealmeans$species=="ACPE")]<-"ACE.PEN"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="ACRU")]<-"ACE.RUB"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="NEMU")]<-"ILE.MUC"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="ILVE")]<-"ILE.VER"
HFrealmeans$GEN.SPA[which(HFrealmeans$species=="VACO")]<-"VAC.COR"
HFrealmeans$scenario<-"field"
HFrealmeans<-dplyr::select(HFrealmeans,-species)
coef(mod.FLS.small)
new.data<-data.frame(GEN.SPA=rep(unique(small$GEN.SPA),9),
                     Force=rep(c(0,1),each=45),
                     Chill=rep(c(.67,1,0),30),
                     Light=rep(c(1),90))



prediction<-predict(mod.FLS.small,newdata=new.data,probs = c(.055,.945))
predy<-cbind(new.data,prediction)


predy$scenario<-NA
predy$scenario[which(predy$Force==0 & predy$Chill==.67)]<-"historic"
predy$scenario[which(predy$Force==1 & predy$Chill==.67)]<-"warm 5"
predy$scenario[which(predy$Force==1 & predy$Chill== 1)]<-"5+chill"
predy$scenario[which(predy$Force==1 & predy$Chill== 0)]<-"5-chill"
unique(predy$scenario)
predy <- na.omit(predy) 
predy<-dplyr::select(predy,GEN.SPA,Estimate,Q5.5,Q94.5,scenario)

predy<-rbind(predy,HFrealmeans)

predy$grouper<-NA
predy$grouper[which(predy$scenario %in% c("field"))]<-"historic-observed"
predy$grouper[which(predy$scenario %in% c("historic"))]<-"historic-predicted"
predy$grouper[which(predy$scenario %in% c("warm 5", "warm 10"))]<-"warm only"
predy$grouper[which(predy$scenario %in% c("5-chill", "10-chill"))]<-"warm,reduce chill"
predy$grouper[which(predy$scenario %in% c("5+chill", "10+chill"))]<-"warm,increase chill"

predy<-filter(predy,GEN.SPA %in% c("ACE.PEN","ACE.RUB","ILE.MUC","ILE.VER","VAC.COR"))
predy %>%
  arrange(Estimate) %>%
  mutate(scenario = factor(scenario, levels=c("field","historic","warm 5","warm 10", "5-chill","10-chill", "5+chill","10+chill"))) %>%
  ggplot(aes(scenario,Estimate))+geom_point(aes(color=grouper))+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0,color=grouper))+facet_wrap(~GEN.SPA)+theme_bw()+geom_hline(yintercept=0)

prey<-filter(predy,scenario %in% c("field","historic"))
prey$scenario[which(prey$scenario %in% c("field"))]<-"observed"
prey$scenario[which(prey$scenario %in% c("historic"))]<-"predicted"

prey$GEN.SPA[which(prey$GEN.SPA=="ACE.PEN")]<-"A. pensylvanicum"
prey$GEN.SPA[which(prey$GEN.SPA=="ACE.RUB")]<-"A. rubrum"
prey$GEN.SPA[which(prey$GEN.SPA=="ILE.MUC")]<-"I. mucronata"
prey$GEN.SPA[which(prey$GEN.SPA=="ILE.VER")]<-"I. veticillata"
prey$GEN.SPA[which(prey$GEN.SPA=="VAC.COR")]<-"V. corymbosum"


jpeg("Plots/fieldmodcomparisions_freescale.jpeg",width = 5, height = 6, units = 'in', res=300)
ggplot(prey,aes(scenario,Estimate))+geom_point(aes())+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0))+facet_wrap(~GEN.SPA,scales="free_y")+theme_bw()+geom_hline(yintercept=0)+ theme(strip.text = element_text(face = "italic"))
dev.off()
jpeg("Plots/fieldmodcomparisions.jpeg",width = 5, height = 6, units = 'in', res=300)
ggplot(prey,aes(scenario,Estimate))+geom_point(aes())+geom_errorbar(aes(ymin=Q5.5,ymax=Q94.5,width=0))+facet_wrap(~GEN.SPA)+theme_bw()+geom_hline(yintercept=0)+ theme(strip.text = element_text(face = "italic"))
dev.off()

