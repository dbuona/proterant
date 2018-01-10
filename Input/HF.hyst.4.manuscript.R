####This is the produces output of harvard forest. Jan 8 2018.
d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
d.spp<-read.csv("hf003-06-mean-spp.csv",header = TRUE)

sum<-group_by(d.spp,species)
sum<-filter(sum,species!=c("QUAL"))
sum<-filter(sum,species!=c("CRSP"))
sum<-filter(sum,species!=c("HAVI"))
sum<-filter(sum,species!=c("TSCA"))
sum<-filter(sum,species!=c("PIST"))
sum<-filter(sum,species!=c("CADE"))

unique(sum$species)
##clean
names(sum)[names(sum) == 'bb.jd'] <- 'leaf budburst'
names(sum)[names(sum) == 'fbb.jd'] <- 'flower budburst'
names(sum)[names(sum) == 'fopn.jd'] <- 'flower open'
names(sum)[names(sum) == 'l75.jd'] <- 'increasing leaf size'

##skip: CADE, CRSP, HAVI,TSCA, pist, qual
sum$species[sum$species=="ACPE"]<-"Acer pensylvanicm"
sum$species[sum$species=="ACRU"]<-"Acer rubrum"
sum$species[sum$species=="ACSA"]<-"Acer saccarum"
sum$species[sum$species=="AMSP"]<-"Amelanchier sp."
sum$species[sum$species=="ARSP"]<-"Aronia sp"
sum$species[sum$species=="BEAL"]<-"Betula allegheniensis"
sum$species[sum$species=="BELE"]<-"Betula lenta"
sum$species[sum$species=="BEPA"]<-"Betula papyrifera"
sum$species[sum$species=="BEPO"]<-"Betula populifolia"
sum$species[sum$species=="COAL"]<-"Cornus alternifolia"
sum$species[sum$species=="FAGR"]<-"Fagus granifolia"
sum$species[sum$species=="FRAM"]<-"Fraxinus americana"
sum$species[sum$species=="ILVE"]<-"Ilex verticillata"
sum$species[sum$species=="KAAN"]<-"Kalmia angustifolia"
sum$species[sum$species=="KALA"]<-"Kalmia latifolia"
sum$species[sum$species=="LYLI"]<-"Lyonia ligustrina"
sum$species[sum$species=="NEMU"]<-"Ilex mucronata"
sum$species[sum$species=="NYSY"]<-"Nyssa sylvatica"
sum$species[sum$species=="COAL"]<-"Cornus alternifolia"
sum$species[sum$species=="POTR"]<-"Populus tremuloides"
sum$species[sum$species=="PRSE"]<-"Prunus seriotina"
sum$species[sum$species=="QURU"]<-"Quercus rubra"
sum$species[sum$species=="QUVE"]<-"Qurcus velutina"
sum$species[sum$species=="RHSP"]<-"Rhododendron sp."
sum$species[sum$species=="SAPU"]<-"Sambucus racemosa"
sum$species[sum$species=="VACO"]<-"Vaccinium corymbosum"
sum$species[sum$species=="VIAL"]<-"Viburnum alnifolium"
sum$species[sum$species=="VICA"]<-"Viburnum cassinoides"

###functional hysteranthy for harvard forest
bb<-gather(sum,phenophase, D.O.Y,3:6)
bb<-filter(bb,phenophase==c("increasing leaf size","flower open"))
meanleaf<-filter(bb,phenophase=="increasing leaf size")
summary(meanleaf)
fun<-ggplot(bb,aes(species,D.O.Y))+stat_summary(aes(shape=phenophase,color=phenophase))+geom_abline(slope=0,intercept=142,color="green")+geom_abline(slope=0,intercept=148,color="dark green")+theme_bw()+ggtitle("Functional definition")+theme(axis.text.x = element_text(angle = 300, hjust = 0))
fun

####physiological hysteranthy
bb<-gather(sum,phenophase,D.O.Y,3:6)
bb2<-filter(bb,phenophase==c("leaf budburst","flower budburst"))
meanlb<-filter(bb,phenophase=="leaf budburst")
summary(meanlb)
phys<-ggplot(bb2,aes(species,D.O.Y))+stat_summary(aes(shape=phenophase,color=phenophase))+geom_abline(slope=0,intercept=119,color="green")+geom_abline(slope=0,intercept=125,color="dark green")+theme_bw()+ggtitle("Physiological definition")+theme(axis.text.x = element_text(angle = 300, hjust = 0))
phys
