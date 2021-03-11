## Started 3 April 2019 ##
## By Cat - based off code from Lizzie's OSPREE plots ##

# Runs from models_stan_plotting.R #

# with COLORS for each species #

muplotfx2 <- function(modelhere,modelhere2,modelhere3, nameforfig, width, height, ylim, xlim, leg1, leg2,leg3,leg4){
  spnum <- unique(dat$GEN.SPA) ### need to change
  setEPS()
  postscript(file.path(figpath, paste("", nameforfig, ".eps", sep="")),
             width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,7,3,10))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main="")
  axis(2, at=1:6, labels=rev(c("Chill", "Photo", "Force","Chill:Photo","Chill:Force","Photo:Force")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  ppeffects <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  for(i in 1:6){#i=1 ## need to change
    pos.y<-(6:1)[i] ### need to change
    pos.x<-(modoutput1[(modoutput1$term==rownameshere[i]),"estimate"]-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"])
    pos.x1<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"])
    pos.x2<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"])
    lines(c(0,pos.x1),rep(pos.y,2),cex=10,lty=1,col="darkblue")
    lines(c(0,pos.x),rep(pos.y-.15,2),cex=10,lty=2,col="darkblue")
    lines(c(0,pos.x2),rep(pos.y-.3,2),cex=10,lty=3,col="darkblue")
    text(0,pos.y, labels="flowering",cex=.1,postion="right")
    text(pos.x,pos.y, labels="budburst",cex=-.1,position="left")
    text(0,pos.y-.15, labels="flowering",cex=-.1,postion="right")
    text(pos.x,pos.y-.15, labels="leafout",cex=-.1,position="left")
    text(0,pos.y-.3, labels="budburst",cex=-.1,position="right")
    text(pos.x,pos.y-.3, labels="leafout",cex=-.1,position="left")
    
  }
}
