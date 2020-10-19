rm(list=ls()) 
options(stringsAsFactors = FALSE)

### makes main effect size plots 
########## THis make a plot#################

library(broom)
library(RColorBrewer)
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
leg.txt <- c("flower","leaf")
source("exp_muplot_brms.R")
source("prep4plot.R")



muplotfx(modelhere,modelhere2, "budburstvsflowering",8, 8, c(0,6), c(-50, 30) , 35, 3.3,35,4)
dev.off()
