####This is the produces output of harvard forest. Jan 8 2018.
d<-read.csv("hf003-05-mean-ind.csv",header=TRUE)
d.spp<-read.csv("hf003-06-mean-spp.csv",header = TRUE)
sum<-group_by(d.spp,species)
names(sum)[names(sum) == 'bb.jd'] <- 'lbb.jd'

###functional hysteranthy for harvard forest
bb<-gather(sum,phenophase,DOY,3:6)
bb<-filter(bb,phenophase==c("l75.jd","fopn.jd"))
meanleaf<-filter(bb,phenophase=="l75.jd")
summary(meanleaf)
fun<-ggplot(bb,aes(species,DOY))+stat_summary(aes(color=phenophase))+geom_abline(slope=0,intercept=142,color="green")+geom_abline(slope=0,intercept=148,color="dark green")+ggtitle("Functional definition")
fun

####physiological hysteranthy
bb<-gather(sum,phenophase,DOY,3:6)
bb2<-filter(bb,phenophase==c("lbb.jd","fbb.jd"))
meanlb<-filter(bb,phenophase=="lbb.jd")
summary(meanlb)
phys<-ggplot(bb2,aes(species,DOY))+stat_summary(aes(color=phenophase))+geom_abline(slope=0,intercept=119,color="green")+geom_abline(slope=0,intercept=125,color="dark green")+ggtitle("Physiological definition")
phys