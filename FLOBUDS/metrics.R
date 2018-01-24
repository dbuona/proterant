##This code tells you how many (and which) species have flowered in each treatment
##To see the most recent metrics just change the filter date in line 9
##started by Dan on 23 Jan 2018
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS")
d<-read.csv("datasheet_alt_format.csv",header=TRUE)
library(tidyverse)
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

