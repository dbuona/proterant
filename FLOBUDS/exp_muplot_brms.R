## Started 3 April 2019 ##
## By Cat - based off code from Lizzie's OSPREE plots ##

# Runs from models_stan_plotting.R #

# with COLORS for each species #

muplotfx <- function(modelhere,modelhere2, nameforfig, width, height, ylim, xlim, leg1, leg2){
  spnum <- unique(dat$GEN.SPA) ### need to change
  pdf(file.path(figpath, paste("", nameforfig, ".pdf", sep="")),
      width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,7,3,10))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main=nameforfig)
  axis(2, at=1:7, labels=rev(c("Intercept","Chill", "Light", "Force","Chill:Light","Chill:Force","Light:Force")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Intercept","b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  ppeffects <- c("b_Intercept","b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  for(i in 1:7){#i=1 ## need to change
    pos.y<-(7:1)[i] ### need to change
    pos.x<-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"]
    lines(modoutput1[(modoutput1$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    pos.x2<-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"]
    lines(modoutput2[(modoutput2$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    points(pos.x,pos.y,cex=1.5,pch=19,col="darkblue")
    points(pos.x2,pos.y,cex=1.5,pch=17,col="darkblue")
    for(spsi in 1:length(spnum)){#
      pos.sps.i<-which(grepl(paste("[",spsi,"]",sep=""),mod.ranef$parameter,fixed=TRUE))
      jitt<-runif(1,0.05,0.5)
      pos.y.sps.i<-pos.y-jitt
      pos.x.sps.i<-mod.ranef[pos.sps.i[i],"mean"]
      pos.x.sps.2<-mod.ranef2[pos.sps.i[i],"mean"]
      lines(mod.ranef[pos.sps.i[i],c("25%","75%")],rep(pos.y.sps.i,2),
            col=alpha(my.pal[spsi], alphahere))
      lines(mod.ranef2[pos.sps.i[i],c("25%","75%")],rep(pos.y.sps.i,2),
            col=alpha(my.pal[spsi], alphahere))
      points(pos.x.sps.i,pos.y.sps.i,cex=0.8,pch=19, col=alpha(my.pal[spsi], alphahere))
      points(pos.x.sps.2,pos.y.sps.i,cex=0.8,pch=17, col=alpha(my.pal[spsi], alphahere))
    }

  }
  par(xpd=TRUE) # so I can plot legend outside
  legend(leg1, leg2, sort(unique(gsub("_", " ", dat$GEN.SPA))), pch=rep(c(19,17),each=10), ### need to change
         col=alpha(my.pal[1:length(spnum)], alphahere),
         cex=1, bty="n", text.font=3)
  
}
