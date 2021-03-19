rm(list=ls()) 
options(stringsAsFactors = FALSE)

### makes main effect size plots 
########## THis make a plot#################
setwd("Documents/git/proterant/FLOBUDS/")
library(broom)
library(RColorBrewer)
library(dplyr)
library(scales)
load("writing/flobud.main.mods.Rda")

figpath <- "Plots"

cols <- adjustcolor("indianred3", alpha.f = 0.3) 

my.pal <- rep(brewer.pal(n = 10, name = "Paired"), 8)
# display.brewer.all()
alphahere = 1

xlab <- "Phenological sensitivity"

spp <- unique(dat$GEN.SPA)



modelhere <-mod.bb.int
modelhere2<-mod.flo.int
modelhere3<-mod.lo.int
leg.txt <- c("flowering","budburst","leafout")
source("exp_muplot_brms.R")
source("prep4plot.R")



muplotfx(modelhere,modelhere2, modelhere3,"all3phases",8, 8, c(0,6), c(-50, 30) , 35, 3.3,35,4)
dev.off()

xlab <- "Difference in phenological sensitivity"
source("weird_muplotsource.R")

png("Plots/grover.png",width = 10, height= 8, units = "in",res = 300)
par(mfrow=c(1,3))
muplotfx2(modoutput1,modoutput2,mod.ranef,mod.ranef2,"Diffplot_bb_flo",8, 8, c(0,6.3), c(-25, 25) , 28, 3.3,35,4,"Flowering","Leaf \nbudburst")


muplotfx3(modoutput3,modoutput2,mod.ranef3,mod.ranef2,"Diffplot_lo_flo",5, 8, c(0,6.3), c(-25, 25) , 28, 3.3,35,4,"Flowering","Leafout")


muplotfx4(modoutput3,modoutput1,mod.ranef3,mod.ranef,"Diffplot_lo_bb",8, 8, c(0,6.3), c(-25, 25) , 27, 3.3,35,4,"Leaf \nBudburst","Leafout")


dev.off()



#### differences


