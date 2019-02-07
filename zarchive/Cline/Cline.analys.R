rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Desktop/Cline")

tree<-read.csv("individual_phenometrics_data.csv",header=TRUE)

table(tree$Common_Name)
library(dplyr)
library(ggplot2)
library(car) 

###Red maple
red<-filter(tree, Species=="rubrum")
table(red$Phenophase_Description)
red<-filter(red, Phenophase_Description==c("Increasing leaf size","Open flowers"))
red<-filter(red, First_Yes_DOY<200)
red<-filter(red, First_Yes_Year>2010)

            

table(red$First_Yes_Year)
red17<-filter(red, First_Yes_Year==2017)
red16<-filter(red, First_Yes_Year==2016)
red15<-filter(red, First_Yes_Year==2015)
red13<-filter(red, First_Yes_Year==2013)


scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red13, main="Red Maple, 2013") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red15, main="Red Maple, 2015") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red16, main="Red Maple, 2016") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red17, main="Red Maple, 2017") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red, main="Red Maple, all years") 
scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=red) 


###red maple with break ing leaf bud
red<-filter(tree, Species=="rubrum")
red<-filter(red, Phenophase_Description==c("Breaking leaf buds","Open flowers"))
red<-filter(red, First_Yes_DOY<200)
red<-filter(red, First_Yes_Year>2010)


table(red$First_Yes_Year)
red17<-filter(red, First_Yes_Year==2017)
red16<-filter(red, First_Yes_Year==2016)
red15<-filter(red, First_Yes_Year==2015)
red13<-filter(red, First_Yes_Year==2013)


scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red13, main="Red Maple, 2013") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red15, main="Red Maple, 2015") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red16, main="Red Maple, 2016") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red17, main="Red Maple, 2017") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red, main="Red Maple, all years") 
scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=red) 

###red maple with break ing leaf bud anf flower bud
red<-filter(tree, Species=="rubrum")
red<-filter(red, Phenophase_Description==c("Breaking leaf buds","Flowers or flower buds"))
red<-filter(red, First_Yes_DOY<200)
red<-filter(red, First_Yes_Year>2010)


table(red$First_Yes_Year)
red17<-filter(red, First_Yes_Year==2017)
red16<-filter(red, First_Yes_Year==2016)
red15<-filter(red, First_Yes_Year==2015)
red13<-filter(red, First_Yes_Year==2013)


scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red13, main="Red Maple, 2013") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red15, main="Red Maple, 2015") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red16, main="Red Maple, 2016") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red17, main="Red Maple, 2017") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=red, main="Red Maple, all years") 

scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=red) 

MA<-filter(red, State=="VA")
scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=MA) 

ß####sugar maple

sug<-filter(tree, Species=="saccharum")
sug<-filter(sug, Phenophase_Description==c("Increasing leaf size","Open flowers"))
sug<-filter(sug, First_Yes_DOY<250)
sug<-filter(sug, Latitude>34)
table(sug$First_Yes_Year)
sug17<-filter(sug, First_Yes_Year==2017)
sug16<-filter(sug, First_Yes_Year==2016)

scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=sug17, main="Sugar Maple 2017") 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=sug16, main="Sugar Maple 2016") 

scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=sug) 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=sug, main="Sugar Maple all years") 
###Yellow birch

bet<-filter(tree, Species=="alleghaniensis")
bet<-filter(bet, Phenophase_Description==c("Increasing leaf size","Open flowers"))
bet<-filter(bet, First_Yes_DOY<250)
bet17<-filter(sug, First_Yes_Year==2017)


scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=bet17, main="Betula 2017") 
scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=bet) 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=bet, main="Betula all years") ###not enough data

#####Cornus
corn<-filter(tree, Species=="florida")
corn<-filter(corn, Phenophase_Description==c("Increasing leaf size","Open flowers"))
corn<-filter(corn, First_Yes_DOY<250)
corn<-filter(corn, Latitude>30)
table(corn$First_Yes_Year)
corn17<-filter(corn, First_Yes_Year==2017)


scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=corn17, main="Dogwood 2017") 
scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=corn) 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=corn, main="Dogwood all years") 


###Red Oak

oak<-filter(tree, Species=="rubra")
oak<-filter(oak, Phenophase_Description==c("Increasing leaf size","Open flowers"))
oak<-filter(oak, First_Yes_DOY<250)
oak<-filter(oak, Latitude>38)
table(oak$First_Yes_Year)
oak17<-filter(oak, First_Yes_Year==2017)


scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=oak17, main="oak 2017") 
scatterplot(First_Yes_DOY~ First_Yes_Year| Phenophase_Description, data=oak) 
scatterplot(First_Yes_DOY~ Latitude| Phenophase_Description, data=oak, main="oak all years") 




