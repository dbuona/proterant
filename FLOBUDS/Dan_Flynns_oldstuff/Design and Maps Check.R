# Creating maps and assigning treatments for Budburst experiment 2016

# This script checks treatment assignments and chamber locations for each cutting.
setwd("~/Documents/git/budchill")

library(gdata) # for read.xls
d <- read.xls("~/Documents/git/budchill/Budburst Chill Datasheet 2016-01-04.xlsx")

allsp <- sort(unique(d[,"sp"]))


# How many individuals are duplicated within new treatments?

table(d$ind, d$New.Treat)

table(d$ind, d$New.Treat.2) # Should all be 1, mistake in randomization.
summary(as.factor(as.vector(table(d$ind, d$New.Treat.2)))) # should have as few 2 and 3 as possible

table(d$ind, d$Final.Treat) # 

table(d$sp, d$Final.Treat) # 

summary(as.factor(as.vector(table(d$ind, d$Final.Treat)))) # at least now fewer than before
hist(as.vector(table(d$ind, d$Final.Treat)))

#labs <- read.csv("~/Documents/git/budchill/design/Avery Labels 3per.csv")
#write.csv(rep(labs$V4, each = 3), file="LabsTreat.csv")

