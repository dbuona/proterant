###provides some additional cleaning for michigan trees dataset. 
mich.tree<-read.tree("pruned_for_mich.tre")
mich.data<-read.csv("mich_data_full.csv")

#clean av fruit time
mich.data$av_fruit_time[mich.data$av_fruit_time=="persistant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="persitant"]<-12
mich.data$av_fruit_time[mich.data$av_fruit_time=="unreported"]<-9    
mich.data$av_fruit_time<-as.numeric(mich.data$av_fruit_time)

mich.data$fruiting<-NA
mich.data$fruiting<-mich.data$av_fruit_time
#mich.data$fruiting[mich.data$fruiting==19]<-7
mich.data$fruiting[mich.data$fruiting=="persistant"]<-12
mich.data$fruiting[mich.data$fruiting=="persitant"]<-12
mich.data$fruiting[mich.data$fruiting=="unreported"]<-9                                      
mich.data$fruiting<-as.numeric(mich.data$fruiting)

mich.data["pro3"]<-NA
mich.data$pro3[mich.data$Phen.sequence == "pro"] <- 1
mich.data$pro3[mich.data$Phen.sequence == "pro/syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "syn/ser"] <- 0
mich.data$pro3[mich.data$Phen.sequence== "ser"] <- 0 
mich.data$pro3[mich.data$Phen.sequence== "hyst"] <- 0

##clean remaining NA's
subset(mich.data, Genus=="Populus")
##coding salix myriocides to be intorlerant of shade as the rest of the genus
mich.data$shade_bin[mich.data$name == "Salix_myricoides"] <- 0
mich.data$fruiting[mich.data$name == "Populus_nigra"] <- 5.5
