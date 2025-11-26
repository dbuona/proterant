###cleaning part II
source("cleaning/cleaning.R")

#############CREATE A DATA SHEET THAT HAS THE THREE FLOWER TYPES IN SEPARATE COLUMNS#########################
#################THIS WILL BE USEFUL FOR COMPARING CHANGES IN PROTANDRY OR PROTOGYNY######################

####find the first day when species reached 15
d.leaf<-filter(dx,leafphase==15)
first<-aggregate(d.leaf$doy.final, by = list(d.leaf$id), min)
####combine with all data
dater<-as.data.frame(unique(dx$id))
colnames(dater)<- c("id")
colnames(first)<-c("id","leaf_day(15)")
dater<-full_join(dater,first,by="id") ###now you have a data set with first leaves



### flowers (currrently mixed only)
unique(dx$flophase)
d.flo<-filter(dx,flophase %in% c("60","60-F","60-M")) ###This is all the 60s ever I have

##added in 2024
d.dich<-filter(d.flo, flophase!="60")
dich<-d.dich %>%group_by(id,flophase) %>%summarise(ff=min(doy.final))

firstflo<-aggregate(d.flo$doy.final, by = list(d.flo$id),min )

colnames(firstflo)<-c("id","flo_day")
dater<-full_join(dater,firstflo,by="id")

L<-filter(dx,leafphase %in% c(9,11,15))
L1<-aggregate(L$doy.final, by = list(L$id), min)
colnames(L1)<-c("id","budburst(9)")
dater<-full_join(dater,L1,by="id")

LLL<-filter(dx,leafphase %in% c(11,15))
L2<-aggregate(LLL$doy.final, by = list(LLL$id), min)
colnames(L2)<-c("id","Leaf_expand(11)")
dater<-full_join(dater,L2,by="id")

treats<-dplyr::select(d,1:7)
treats<-distinct(treats)

test<-left_join(dater,treats)
write.csv(test,"flobudsdata.use.csv",row.names = FALSE)
