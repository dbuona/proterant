###randomizing code
dat<-read.csv("pre_sample_datasheet.csv") ###this is the data

##assign divide the species between treatments
Z <- block_ra(block_var = dat$species, condition_names = c("WL0", "WS0", "WL1","WS1","CL0","CS0","CL1","CS1"))

###add assignemnt to spreadsheet
dat$assignment<-Z

###check it###############################
Ac<-filter(dat,species=="rubrum")
table(Ac$assignment)
cp<-filter(dat,species=="peregrina")
table(cp$assignment)
cc<-filter(dat,species=="cornuta")
table(cc$assignment)
### It worked

###assign to beakers
X <- block_ra(block_var = dat$assignment,condition_names = 1:24)
dat$group<-X
table(dat$group)

##did the properly assign species withing treatments to groups of 3
WS0<-filter(dat,assignment=="WS0")
table(WS0$group) ###yes but there are extras

