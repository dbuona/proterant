### make circle trees for final manuscript resubmission
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/proterant/Input")


library(ggtree,quietly=TRUE)
library(gridExtra,quietly=TRUE)

mich.data<-read.csv("datasheets_derived/MTSV_USFS/michdata_final.csv")
mich.tre<-read.tree("datasheets_derived/MTSV_USFS/michtre_final.tre")

silv.data<-read.csv("datasheets_derived/MTSV_USFS/silvdata_final.csv")
silv.tre<-read.tree("datasheets_derived/MTSV_USFS/silvtre_final.tre")


mich.data$synth<-mich.data$pro2+mich.data$pro3
mich.data$synth2[which(mich.data$synth==0)] <- "always seranthous"
mich.data$synth2[which(mich.data$synth==1)] <- "transitional"
mich.data$synth2[which(mich.data$synth==2)] <- "always hysteranthous"

tr <-mich.tre

par(mar=c(0,0,0,0))

dd <- data.frame(taxa  = mich.data$name, hysteranthy = mich.data$synth2 )

p<-ggtree(tr,layout="circular")
p <- p %<+% dd +geom_tippoint(aes(color=hysteranthy,shape=hysteranthy),size=2)
p<-p+theme(legend.position = "none")



silv.data$synth<-silv.data$pro2+silv.data$pro3
silv.data$synth2[which(silv.data$synth==0)] <- "always seranthous"
silv.data$synth2[which(silv.data$synth==1)] <- "transitional"
silv.data$synth2[which(silv.data$synth==2)] <- "always hysteranthous"

#### make a tree for USFS svics

tr1 <-silv.tre

dd1 <- data.frame(taxa  = silv.data$name, hysteranthy = silv.data$synth2 )


q<-ggtree(tr1,layout="circular")
q <- q %<+% dd1 +geom_tippoint(aes(color=hysteranthy,shape=hysteranthy),size=2)

q<-q+theme(legend.position="bottom")


setEPS()
postscript("..//sub_projs/cicletrees.eps",width = 8, height = 4)
ggpubr::ggarrange(p,q,ncol=2,legend="bottom",common.legend = TRUE,labels = c("a)","b)"))
dev.off()
                  