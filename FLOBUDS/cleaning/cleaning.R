###This is the main cleaning file for flobuds data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/FLOBUDS")
d<-read.csv("datasheet_alt_format.csv",header=TRUE)
library(tidyverse)
library(lubridate)
library("Hmisc")

###Clean date
unique(d$date)
d$date[d$date=="26-Nov"]<-"11/26/17"
d$Date<-d$date
d<-separate(d,date,c("month","day","year"))

d$Date<-as.Date(d$Date,format =  "%m/%d/%y")

d$doy<-yday(d$Date)
unique(d$doy)
start<-yday("2017/11/21")
365-325
d$doy.adjusted<-ifelse(d$year==17,d$doy-start,40+(d$doy))
unique(d$doy.adjusted)

#cold treatment wen in Dec/19/2017
start2<-yday("2017/12/19")
365-353
d$doy.adjusted2<-ifelse(d$year==17,d$doy-start2,12+(d$doy))
unique(d$doy.adjusted2)

d$doy.final<-ifelse(d$Chill==0,d$doy.final<-d$doy.adjusted,d$doy.final<-d$doy.adjusted2)

d<-filter(d,doy.final>=0)
unique(d$doy.final)
max(d$doy.final)
q<-filter(d,Chill==1)
r<-filter(d,Chill==0)
unique(q$doy.final)
unique(r$doy.final)
###give each entry a unique id
d<-unite(d,id,name,good_flaskid,sep="_",remove = FALSE)

###clean treatment
d<-unite(d,treatcode,Force,Light,Chill,sep="_",remove = FALSE)
d$Light<-ifelse(d$treatcode=="W_S_0","L",d$Light)
d$Light<-ifelse(d$treatcode=="W_L_0","S",d$Light)
d<-dplyr::select(d,-treatcode)## we'll add this later, but in the meantime it will mess up all the gather() commands if we dont drop it

###clean species
colnames(d)
p<-filter(d, GEN.SPA=="PRU.PEN")
unique(p$name)
d$GEN.SPA[d$name=="PRUVIR 7 HF DB"]<-"AME.SPP"
d$GEN.SPA[d$name=="PRUPEN 21 HF DB"]<-"BET.SPP"


#### gather locations so it is informative but not catagorical
dx<-gather(d,"flower_location","flophase",c(8,10))
dx<-gather(dx,"leaf_location","leafphase",c(8,9))


####clean extra data points from first experiemnt
dx<-filter(dx,doy.final<=112)

####clean flowering
unique(dx$flophase)####this needs alot of cleaning!!!!
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
###more cleaning
dx$flophase[dx$flophase=="0,69-F,67-F"]<-"0,69-F,67-M"
dx$flophase[dx$flophase=="0,065-F,0"]<-"0,65-F,0"
dx$flophase[dx$flophase=="0,60-F,65-N"]<-"0,60-F,65-M"
dx$flophase[dx$flophase=="0,60-F,060-M"]<-"0,60-F,60-M"
dx$flophase[dx$flophase=="0,69-F,690M"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="0,69-M,69-M"]<-"0,69-F,69-M"
dx$flophase[dx$flophase=="0,60-F,60M"]<-"0,60-F,60-M"
dx$flophase[dx$flophase=="0,69-F,065-M"]<-"0,69-F,65-M"
unique(dx$flophase)


dx<-separate(dx,flophase,c("mixphase","femphase","malephase"),sep=",")
dx<-gather(dx,flotype,flophase,17:19)
unique(dx$flophase)


###### make everything 60######## for computation sake this way if first flower was score at 65 (etc) it makes the analysis.
dx$flophase[dx$flophase=="61"]<-60
dx$flophase[dx$flophase=="62"]<-60
dx$flophase[dx$flophase=="63"]<-60
dx$flophase[dx$flophase=="64"]<-60
dx$flophase[dx$flophase=="65"]<-60
dx$flophase[dx$flophase=="67"]<-60

dx$flophase[dx$flophase=="60-F"]<-"60-F"
dx$flophase[dx$flophase=="61-F"]<-"60-F"
dx$flophase[dx$flophase=="62-F"]<-"60-F"
dx$flophase[dx$flophase=="63-F"]<-"60-F"
dx$flophase[dx$flophase=="64-F"]<-"60-F"
dx$flophase[dx$flophase=="65-F"]<-"60-F"
dx$flophase[dx$flophase=="67-F"]<-"60-F"

dx$flophase[dx$flophase=="60-M"]<-"60-M"
dx$flophase[dx$flophase=="61-M"]<-"60-M"
dx$flophase[dx$flophase=="62-M"]<-"60-M"
dx$flophase[dx$flophase=="63-M"]<-"60-M"
dx$flophase[dx$flophase=="64-M"]<-"60-M"
dx$flophase[dx$flophase=="65-M"]<-"60-M"
dx$flophase[dx$flophase=="67-M"]<-"60-M"
dx$flophase[dx$flophase=="60M"]<-"60-M"
table(dx$flophase)

### clean flophase
dx$flotype<-ifelse(dx$GEN.SPA=="VAC.COR","mixphase",dx$flotype)
dx$flotype<-ifelse(dx$GEN.SPA=="PRUN.PEN","mixphase",dx$flotype)
dx$flotype<-ifelse(dx$GEN.SPA=="PRU.VIR","mixphase",dx$flotype)

##Comper
dx$flophase<-ifelse(dx$id=="COMPER 1 HF DB_CL1_24" &dx$flophase==60,"60-F",dx$flophase)
dx$flophase<-ifelse(dx$id=="COMPER 6 HF DB_WS0_21" &dx$flophase==60,"60-F",dx$flophase)
dx$flophase<-ifelse(dx$id=="COMPER 2 HF DB_CL0_1" &dx$flophase==60,"60-F",dx$flophase)


#to do:
#put in treatment values (IE 8 and 12 photoperiod)?
#Calculate chilling?
#restrict to 112 days for both chilling treatment

