####Duration of flowering
#### overlap of flowering and leaftime
3###May 9 2018

setwd("~/Documents/git/proterant/FLOBUDS")
source("cleaning.R")
dater<-as.data.frame(unique(dx$id))
colnames(dater)<- c("id")
unique(dx$flophase)

d.flo<-filter(dx,flophase %in% c(60))

firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),min )
lastflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),max )
colnames(firstflo)<-c("id","flo_start")
colnames(lastflo)<-c("id","flo_end")
dur.dat<-full_join(firstflo,lastflo)
dur.dat<-full_join(dater,dur.dat)
dur.dat$duration<-dur.dat$flo_end-dur.dat$flo_start
dur.dat<-na.omit(dur.dat)
dur.dat$sex="bisexual"
d.flo<-filter(dx,flophase %in% c("60-F"))

firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),min )
lastflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),max )
colnames(firstflo)<-c("id","flo_start")
colnames(lastflo)<-c("id","flo_end")
dur.dat1<-full_join(firstflo,lastflo)
dur.dat1<-full_join(dater,dur.dat1)
dur.dat1$duration<-dur.dat1$flo_end-dur.dat1$flo_start
dur.dat1<-na.omit(dur.dat1)
dur.dat1$sex<-"female"
D<-rbind(dur.dat,dur.dat1)

d.flo<-filter(dx,flophase %in% c("60-M"))

firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),min )
lastflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),max )
colnames(firstflo)<-c("id","flo_start")
colnames(lastflo)<-c("id","flo_end")
dur.dat2<-full_join(firstflo,lastflo)
dur.dat2<-full_join(dater,dur.dat2)
dur.dat2$duration<-dur.dat2$flo_end-dur.dat2$flo_start
dur.dat2<-na.omit(dur.dat2)
dur.dat2$sex<-"male"
DATA<-rbind(D,dur.dat2)

dxx<-dplyr::select(dx,1:7)
dxx<-distinct(dxx)
goober<-left_join(DATA,dxx, by="id")
goober<-unite(goober, treatment, Chill, Force, Light, remove=FALSE)

ggplot(goober, aes(treatment,duration))+stat_summary(aes(color=sex))+facet_wrap(~GEN.SPA)

       