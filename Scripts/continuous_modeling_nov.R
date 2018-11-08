
###individuals with more than 10 observations
checker1<-d %>% group_by(tree.id) %>% summarise(non_na_count = sum(!is.na(offset.funct)))
checker1<-filter(checker1, non_na_count>=10)
table(checker1$tree.id)
use.id<-checker1$tree.id
d.plus<-filter(d, tree.id %in% use.id)

##species to use
###ACRU, QURU, ACPE, most complete observations, each hysteranthy class

###species with more than 45 total boservations
use.sp<-c("ACPE","ACRU","QURU")
#"KALA","BEAL","AMSP","VACO

d.plus<-filter(d.plus, species %in% use.sp)

ggplot(d.plus)+stat_summary(fun.data = "mean_cl_boot",geom="errorbar",aes(tree.id, offset.phys,color=species),linetype="solid")+stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",aes(tree.id,offset.funct,color=species),linetype="dashed")+stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",aes(tree.id,offset.inter,color=species),linetype="dotdash")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))

m<-lmer(offset.inter~tree.id+(1|species)+(1|year),data=d.plus)
n<-lmer(offset.phys~tree.id+(1|species)+(1|year),data=d.plus)
o<-lmer(offset.funct~tree.id+(1|species)+(1|year),data=d.plus)
Anova(m)
Anova(n)
Anova(o)
####This indicates that idividuals vary in the phys and inter, but bot
#### interannual


ggplot(d.plus)+geom_line(aes(year,offset.funct, group=tree.id, color=species),linetype="dashed")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+geom_line(aes(year,offset.phys, group=tree.id, color=species),linetype="solid")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))+geom_line(aes(year,offset.inter, group=tree.id, color=species),linetype="dotdash")+theme_bw()+scale_color_manual(values=c("darkgreen","red","blue"))

mm<-lmer(offset.inter~year+(1|tree.id),data=d.plus)
nn<-lmer(offset.phys~year+(1|tree.id),data=d.plus)
oo<-lmer(offset.funct~year+(1|tree.id),data=d.plus)
summary(mm)
Anova(mm)


  ###This is the continuous model not using michigan tree so we can address other hypothesese







setdiff(d$species, traits$species)

dater<-left_join(d,traits, by="species")
dater$cent_flo_time<-(dater$fopn.jd-mean(dater$fopn.jd,na.rm=TRUE))/(2*sd(dater$fopn.jd,na.rm=TRUE))

##you should take a mean flowering time for each spcecies

ggplot(dater, aes(offset))+geom_histogram()

goob<-lmer(offset~cent_pol+seed_cent+(1|tree.id)+(1|year),data=dater)
summary(goob)
library(brms)
brm(offset~pol+(1|year), data=dater, family=gaussian())
