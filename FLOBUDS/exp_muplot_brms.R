## Started 3 April 2019 ##
## By Cat - based off code from Lizzie's OSPREE plots ##

# Runs from models_stan_plotting.R #

# with COLORS for each species #

muplotfx <- function(modelhere,modelhere2, nameforfig, width, height, ylim, xlim, leg1, leg2,leg3,leg4){
  spnum <- unique(dat$GEN.SPA) ### need to change
  setEPS()
  postscript(file.path(figpath, paste("", nameforfig, ".eps", sep="")),
      width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,7,3,10))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main="")
  axis(2, at=1:6, labels=rev(c("Chill", "Light", "Force","Chill:Light","Chill:Force","Light:Force")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  ppeffects <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  for(i in 1:6){#i=1 ## need to change
    pos.y<-(6:1)[i] ### need to change
    pos.x<-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"]
    lines(modoutput1[(modoutput1$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    pos.x2<-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"]
    lines(modoutput2[(modoutput2$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    points(pos.x,pos.y,cex=1.5,pch=17,col="darkblue")
    points(pos.x2,pos.y,cex=1.5,pch=19,col="darkblue")
    for(spsi in 1:length(spnum)){#
      pos.sps.i<-which(grepl(paste("[",spsi,"]",sep=""),mod.ranef$parameter,fixed=TRUE))
      jitt<-0.075
      pos.y.sps.i<-pos.y-jitt*c(spnum[2])
      print(pos.y.sps.i)
      pos.x.sps.i<-mod.ranef[pos.sps.i[i],"mean"]
      pos.x.sps.2<-mod.ranef2[pos.sps.i[i],"mean"]
      lines(mod.ranef[pos.sps.i[i],c("25%","75%")],rep(pos.y.sps.i,2),
            col=alpha(my.pal[spsi], alphahere))
      lines(mod.ranef2[pos.sps.i[i],c("25%","75%")],rep(pos.y.sps.i,2),
            col=alpha(my.pal[spsi], alphahere))
      points(pos.x.sps.i,pos.y.sps.i,cex=0.8,pch=17, col=alpha(my.pal[spsi], alphahere))
      points(pos.x.sps.2,pos.y.sps.i,cex=0.8,pch=19, col=alpha(my.pal[spsi], alphahere))
    }

  }
  par(xpd=TRUE) # so I can plot legend outside
  legend(leg1, leg2, sort(unique(gsub("_", " ", c("Comptonia perigrina","Viburnum acerifolium", "Ilex verticillata","Vaccinium corymbosum","Acer penyslyvanicum",
  "Prunus virginiana","Corylus cornuta","Prunus pensylvanica","Ilex mucronata","Acer rubrum")))), pch=rep(c(15),each=10), ### need to change
         col=alpha(my.pal[1:length(spnum)], alphahere),
         cex=.8, bty="n", text.font=3)
  legend(leg3, leg4, legend=leg.txt, pch=c(19,17), ### need to change
         cex=.8, bty="n", text.font=3)
  
}
