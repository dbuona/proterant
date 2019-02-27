
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS")
source("cleaning/cleaning.R")

d.leaf<-filter(dx,leafphase %in% c(7,9))
d.leafy<- d.leaf %>% group_by(id) %>% summarise(LBB=min(doy.final))



d.leaf3<-filter(dx,leafphase %in% c(11))
d.leafy3<- d.leaf3 %>% group_by(id) %>% summarise(L11=min(doy.final))

d.leaf4<-filter(dx,leafphase %in% c(15))
d.leafy4<- d.leaf4 %>% group_by(id) %>% summarise(L15=min(doy.final))


dater<-as.data.frame(unique(dx$id))
colnames(dater)<- c("id")

dater<-left_join(dater,d.leafy)
dater<-left_join(dater,d.leafy3)
dater<-left_join(dater,d.leafy4)

##add flowering
d.flo<-filter(dx,flophase %in% c("60","60-F","60-M")) ###This is all the 60s ever I have

### get rid of sex specific descriptor because that is alread in a colume
d.flo$flophase<-"60"
d.floy<- d.flo %>% group_by(id) %>% summarise(F60=min(doy.final))
dater<-left_join(dater,d.floy)

dxx<-dplyr::select(dx,1:7)
colnames(dxx)
dxx<-distinct(dxx)
dater<-left_join(dxx,dater)
write.csv(dater,"Flo_buds_for_thermoperiodicity.csv",row.names=FALSE)

