###code for plotting sims
## Started 3 April 2019 ##
## By Cat - based off code from Lizzie's OSPREE plots ##

# Runs from models_stan_plotting.R #

# with COLORS for each species #

muplotfx <- function(modelhere,modelhere2, nameforfig, width, height, ylim, xlim){
 ### need to change
  setEPS()
  postscript(file.path(figpath, paste("", nameforfig, ".eps", sep="")),
      width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,7,3,10))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main=nameforfig)
  axis(2, at=1:3, labels=rev(c("Intercept","Force","Chill")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Intercept","b_temp", "b_chill") ### Need to change
  ppeffects <- c("b_Intercept","b_temp", "b_chill") ### Need to change
  for(i in 1:3){#i=1 ## need to change
    pos.y<-(3:1)[i] ### need to change
    pos.x<-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"]
    lines(modoutput1[(modoutput1$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    pos.x2<-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"]
    lines(modoutput2[(modoutput2$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    points(pos.x,pos.y,cex=1,pch=19,col="black")
    points(pos.x2,pos.y,cex=1,pch=17,col="black")
    
  }}

 