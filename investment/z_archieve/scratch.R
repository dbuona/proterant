
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

setwd("~/Documents/git/proterant")
info<-read.csv("investment/herbaria_prunus_rec.csv")
colnames(info)
shortlist<-dplyr::filter(info,county!="")
table(shortlist$specificEpithet)
table(info$specificEpithet)

unique(shortlist$specificEpithet)
unique(info$specificEpithet)

testcase<-dplyr::filter(info,specificEpithet %in% c("mexicana","ilicifolia"))
table(testcase$specificEpithet)
write.csv(testcase,"investment/twospecies.csv")
