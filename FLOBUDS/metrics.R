##This code tells you how many (and which) species have flowered in each treatment
##To see the most recent metrics just change the filter date in line 9
##started by Dan on 23 Jan 2018

###Dan added cleaning to the bottom of this.

rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS")
d<-read.csv("datasheet_alt_format.csv",header=TRUE)
library(tidyverse)
library(lubridate)
today<-filter(d,date=="1/23/18")
######CSO###############
CS0<- subset(today, today$Force=="C" & today$Chill=="0" & today$Light=="S")
goo<-subset(CS0, CS0$term.F!=0 | CS0$lat.F!=0)
goo2<-subset(CS0, CS0$term.F!=0 | CS0$lat.F!=0 |CS0$term.L!=0 | CS0$lat.L!=0)
goo3<-subset(CS0, CS0$term.F==0 & CS0$lat.F==0 &CS0$term.L==0 & CS0$lat.L==0)
nrow(goo) #4
nrow(goo2) #10
nrow(goo3) #59


###CL0###########
CL0<- subset(today, today$Force=="C" & today$Chill=="0" & today$Light=="L")
doo<-subset(CL0, CL0$term.F!=0 | CL0$lat.F!=0)
doo2<-subset(CL0, CL0$term.F!=0 | CL0$lat.F!=0 |CL0$term.L!=0 | CL0$lat.L!=0)
doo3<-subset(CL0, CL0$term.F==0 & CL0$lat.F==0 &CL0$term.L==0 & CL0$lat.L==0)

nrow(doo) #9
table(doo$GEN.SPA)
nrow(doo2) #16
table(doo2$GEN.SPA)
nrow(doo3) #56

#####WLO ps this is really "short" daylight treatment

WL0<- subset(today, today$Force=="W" & today$Chill=="0" & today$Light=="L")
qoo<-subset(WL0, WL0$term.F!=0 | WL0$lat.F!=0)
qoo2<-subset(WL0, WL0$term.F!=0 | WL0$lat.F!=0 |WL0$term.L!=0 | WL0$lat.L!=0)
qoo3<-subset(WL0, WL0$term.F==0 & WL0$lat.F==0 &WL0$term.L==0 & WL0$lat.L==0)

nrow(qoo)#9
table(qoo$GEN.SPA)
nrow(qoo2) #31
table(qoo2$GEN.SPA)
nrow(doo3) 

###WS0
WS0<- subset(today, today$Force=="W" & today$Chill=="0" & today$Light=="S")
boo<-subset(WS0, WS0$term.F!=0 | WS0$lat.F!=0)
boo2<-subset(WS0, WL0$term.F!=0 | WS0$lat.F!=0 |WS0$term.L!=0 | WS0$lat.L!=0)
boo3<-subset(WS0, WS0$term.F==0 & WS0$lat.F==0 &WS0$term.L==0 & WS0$lat.L==0)

nrow(boo) #17
table(boo$GEN.SPA)
nrow(boo2) #44
table(boo2$GEN.SPA)
nrow(boo3) 


#####WLO ps this is really "short" daylight treatment

WL0<- subset(today, today$Force=="W" & today$Chill=="0" & today$Light=="L")
qoo<-subset(WL0, WL0$term.F!=0 | WL0$lat.F!=0)
qoo2<-subset(WL0, WL0$term.F!=0 | WL0$lat.F!=0 |WL0$term.L!=0 | WL0$lat.L!=0)
qoo3<-subset(WL0, WL0$term.F==0 & WL0$lat.F==0 &WL0$term.L==0 & WL0$lat.L==0)

##################MORE CHILLING

#####WL1
WL0<- subset(today, today$Force=="W" & today$Chill=="1" & today$Light=="L")
qoo<-subset(WL0, WL0$term.F!=0 | WL0$lat.F!=0)
qoo2<-subset(WL0, WL0$term.F!=0 | WL0$lat.F!=0 |WL0$term.L!=0 | WL0$lat.L!=0)
qoo3<-subset(WL0, WL0$term.F==0 & WL0$lat.F==0 &WL0$term.L==0 & WL0$lat.L==0)

nrow(qoo)#21
table(qoo$GEN.SPA)
nrow(qoo2) #56
table(qoo2$GEN.SPA)
nrow(qoo3) 

###WS1
WS0<- subset(today, today$Force=="W" & today$Chill=="1" & today$Light=="S")
boo<-subset(WS0, WS0$term.F!=0 | WS0$lat.F!=0)
boo2<-subset(WS0, WS0$term.F!=0 | WS0$lat.F!=0 |WS0$term.L!=0 | WS0$lat.L!=0)
boo3<-subset(WS0, WS0$term.F==0 & WS0$lat.F==0 &WS0$term.L==0 & WS0$lat.L==0)

nrow(boo) #20
table(boo$GEN.SPA)
nrow(boo2) #40
table(boo2$GEN.SPA)
nrow(boo3) 

###CL1###
CL0<- subset(today, today$Force=="C" & today$Chill=="1" & today$Light=="L")
doo<-subset(CL0, CL0$term.F!=0 | CL0$lat.F!=0)
doo2<-subset(CL0, CL0$term.F!=0 | CL0$lat.F!=0 |CL0$term.L!=0 | CL0$lat.L!=0)
doo3<-subset(CL0, CL0$term.F==0 & CL0$lat.F==0 &CL0$term.L==0 & CL0$lat.L==0)

nrow(doo) #16
table(doo$GEN.SPA)
nrow(doo2) #44
table(doo2$GEN.SPA)
nrow(doo3) #28

CS0<- subset(today, today$Force=="C" & today$Chill=="1" & today$Light=="S")
goo<-subset(CS0, CS0$term.F!=0 | CS0$lat.F!=0)
goo2<-subset(CS0, CS0$term.F!=0 | CS0$lat.F!=0 |CS0$term.L!=0 | CS0$lat.L!=0)
goo3<-subset(CS0, CS0$term.F==0 & CS0$lat.F==0 &CS0$term.L==0 & CS0$lat.L==0)
nrow(goo) #11
table(goo$GEN.SPA)
nrow(goo2) #27
nrow(goo3) #59

### Current flowering 
#W"S"0=17
#W"L"0=9
#CL0=9
#CS0=4
#WL1=21
#WS1=20
#CL1=16
#CS0=11

####all flowering so far
vinners<-subset(today, today$term.F!=0 | today$lat.F!=0)
table(vinners$GEN.SPA)
leafers<-subset(today, today$term.L!=0 | today$lat.L!=0)
nrow(leafers)
table(leafers$GEN.SPA)
##########################################################
##########Take your data, make it short form

d<-read.csv("datasheet_alt_format.csv",header=TRUE)
###Clean date
unique(d$date)
d$date[d$date=="26-Nov"]<-"11/26/17"
d$Date<-d$date
d<-separate(d,date,c("month","day","year"))

d$Date<-as.Date(d$Date,format =  "%m/%d/%y")

d$doy<-yday(d$Date)
unique(d$doy)
start<-yday("2017/11/21")
d$doy.adjusted<-ifelse(d$year==17,d$doy-start,40+(d$doy))
unique(d$doy.adjusted)

###give each entry a unique id
d<-unite(d,id,name,good_flaskid,sep="_",remove = FALSE)

#### gather locations so it is informative but not catagorical
dx<-gather(d,"flower_location","flophase",c(8,10))
dx<-gather(dx,"leaf_location","leafphase",c(8,9))

unique(dx$flophase)####this needs alot of cleaning,but ignore for now
###79 flowering catagories to clean
###change all to format: mixed,female,male
dx$flophase[dx$flophase=="0"]<-"0,0,0"
dx$flophase[dx$flophase=="-"]<-"0,0,0"
dx$flophase[dx$flophase=="9"]<-"9,0,0"
dx$flophase[dx$flophase=="7"]<-"7,0,0"
dx$flophase[dx$flophase=="1"]<-"1,0,0"
dx$flophase[dx$flophase=="NA"]<-"NA,NA,NA"
dx$flophase[is.na(dx$flophase)]<-"NA,NA,NA"
dx$flophase[dx$flophase==""]<-"0,0,0"
dx$flophase[dx$flophase=="na"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="n"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="55"]<-"55,0,0"
dx$flophase[dx$flophase=="51"]<-"51,0,0"
dx$flophase[dx$flophase=="60"]<-"60,0,0"
dx$flophase[dx$flophase=="61"]<-"61,0,0"
dx$flophase[dx$flophase=="63"]<-"63,0,0"
dx$flophase[dx$flophase=="64"]<-"64,0,0"
dx$flophase[dx$flophase=="62"]<-"62,0,0"
dx$flophase[dx$flophase=="65"]<-"65,0,0"
dx$flophase[dx$flophase=="60-M"]<-"0,0,60-M"
dx$flophase[dx$flophase=="67"]<-"0,0,67"
dx$flophase[dx$flophase=="1-M"]<-"0,0,1-M"
dx$flophase[dx$flophase=="51-M"]<-"0,0,51-M"
dx$flophase[dx$flophase=="X"]<-"X,X,X"
dx$flophase[dx$flophase=="55-M"]<-"0,0,55-M"
dx$flophase[dx$flophase=="67-M"]<-"0,0,67-M"
dx$flophase[dx$flophase=="69"]<-"69,0,0"
dx$flophase[dx$flophase=="62-M"]<-"0,0,62-M"
dx$flophase[dx$flophase=="65-M"]<-"0,0,65-M"
dx$flophase[dx$flophase=="65-M,60-F"]<-"0,60-F,65-M"
dx$flophase[dx$flophase=="67,9"]<-"67,0,0"
dx$flophase[dx$flophase=="69,51"]<-"69,0,0"
dx$flophase[dx$flophase=="65-M,65-F"]<-"0,65-F,65-M"
dx$flophase[dx$flophase=="69,9"]<-"69,0,0"
dx$flophase[dx$flophase=="69,55"]<-"69,0,0"
dx$flophase[dx$flophase=="67-M,F"]<-"0,67-F,67-M"
dx$flophase[dx$flophase=="69-M"]<-"0,0,69-M"
dx$flophase[dx$flophase=="667-M"]<-"0,0,67-M"
dx$flophase[dx$flophase=="11"]<-"11,0,0"
dx$flophase[dx$flophase=="D"]<-"D,D,D"
dx$flophase[dx$flophase=="63-M"]<-"0,0,63-M"
dx$flophase[dx$flophase=="60-F"]<-"0,60-F,0"
dx$flophase[dx$flophase=="60-F"]<-"0,60-F,0"
dx$flophase[dx$flophase=="Na"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="69-"]<-"69?,69?,?"
dx$flophase[dx$flophase=="lost"]<-"0,0,lost"
dx$flophase[dx$flophase=="60-m"]<-"0,0,60-M"
dx$flophase[dx$flophase=="?"]<-"?,?,?"
dx$flophase[dx$flophase=="61-M"]<-"0,0,61-M"
dx$flophase[dx$flophase=="60-F"]<-"0,59-F,0"
dx$flophase[dx$flophase=="55-F"]<-"0,55-F,0"
dx$flophase[dx$flophase=="61-F"]<-"0,61-F,0"
dx$flophase[dx$flophase=="63-F"]<-"0,63-F,0"
dx$flophase[dx$flophase=="65-F"]<-"0,65-F,0"
dx$flophase[dx$flophase=="65%"]<-"65,0,0"
dx$flophase[dx$flophase=="51-F"]<-"0,51-F,0"
dx$flophase[dx$flophase=="62-F"]<-"0,62-F,0"
dx$flophase[dx$flophase=="51?"]<-"51?,0,0"
dx$flophase[dx$flophase=="67-F"]<-"0,67-F,0"
dx$flophase[dx$flophase=="65-F,55-M"]<-"0,65-F,55-M"
dx$flophase[dx$flophase=="N"]<-"NA,NA,NA"
dx$flophase[dx$flophase=="67-F,63-M"]<-"0,67-F,63-M"
dx$flophase[dx$flophase=="67-F,65-M"]<-"0,67-F,65-M"
dx$flophase[dx$flophase=="69-F"]<-"0,69-F,0"
dx$flophase[dx$flophase=="69-F,51-M"]<-"0,69-F,51-M"
dx$flophase[dx$flophase=="69-F,67-M"]<-"0,69-F,67-M"
dx$flophase[dx$flophase=="65-F,51-M"]<-"0,65-F,51-M"
dx$flophase[dx$flophase=="69-F,51"]<-"0,69-F,51-M"
dx$flophase[dx$flophase=="69-F,67"]<-"0,69-F,67-M"
dx$flophase[dx$flophase=="69-M,60-F"]<-"0,60-F,69-M"
dx$flophase[dx$flophase=="65-m"]<-"0,0,65-M"
dx$flophase[dx$flophase=="69-F,60-M"]<-"0,69-F,60-M"
dx$flophase[dx$flophase=="69-F,69-M"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="60-M,69-F"]<-"0,69-F,60-M"
dx$flophase[dx$flophase=="69-M,65-F"]<-"0,65-F,69-M"
dx$flophase[dx$flophase=="67-F,60-M"]<-"0,67-F,60-M"
dx$flophase[dx$flophase=="60-F,55-M"]<-"0,60-F,55-M"
dx$flophase[dx$flophase=="69-M,69-F"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="69-M,67-F"]<-"0,67-F,69-M"
dx$flophase[dx$flophase=="69-M,57-M"]<-"0,67-F,69-M"
dx$flophase[dx$flophase=="69-F,65-M"]<-"0,69-F,65-M"
dx$flophase[dx$flophase=="69-F,65-M"]<-"0,69-F,65-M"
dx$flophase[dx$flophase=="15"]<-"15?,0,0"
dx$flophase[dx$flophase=="59-F"]<-"0,59-F,0"
unique(dx$flophase)

dx<-separate(dx,flophase,c("mixphase","femphase","malephase"),sep=",")
dx<-gather(dx,flotype,flophase,15:17)

######find the first day when species reached 15
d.leaf<-filter(dx,leafphase==15)
first<-aggregate(d.leaf$doy.adjusted, by = list(d.leaf$id), min)

####combine with all data

dater<-as.data.frame(unique(dx$id))
colnames(dater)<- c("id")
colnames(first)<-c("id","leaf_day")
dater<-full_join(dater,first,by="id") ###now you have a data set with first leaves

### flowers (currrently mixed only)
d.flo<-filter(dx,flophase==60)
firstflo<-aggregate(d.flo$doy.adjusted, by = list(d.flo$id), min)
colnames(firstflo)<-c("id","flo_day")
dater<-full_join(dater,firstflo,by="id")

###add female
d.flo<-filter(dx,flophase=="60-F")
firstflo<-aggregate(d.flo$doy.adjusted, by = list(d.flo$id), min)
colnames(firstflo)<-c("id","flo_dayF")
dater<-full_join(dater,firstflo,by="id")

###add male
d.flo<-filter(dx,flophase=="60-M")
firstflo<-aggregate(d.flo$doy.adjusted, by = list(d.flo$id), min)
colnames(firstflo)<-c("id","flo_dayM")
dater<-full_join(dater,firstflo,by="id")
