muplotfx_field <- function(modelhere, nameforfig, width, height, ylim, xlim, leg1, leg2){
  spnum <- unique(daters.fullrecords) ### need to change
  pdf(file.path(figpath, paste("", nameforfig, ".pdf", sep="")),
      width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,7,3,10))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main=nameforfig)
  axis(2, at=1:7, labels=rev(c("Intercept","1990","1991","1992","1993","1994","1995","1996")),las=1)#,"1997","1998","1999","2000","2001","2002","2003","2004","2005"),), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Intercept","b_End_year1990","b_End_year1991", "b_End_year1992", "b_End_year1993","b_End_year1994","b_End_year1995","b_End_year1996")#,"b_1997","b_1998","b_1999","b_2000","b_2001","b_2002","b_2003","b_2004","b_2005") ### Need to change
  ppeffects <- c("b_Intercept","b_End_year1990","b_End_year1991", "b_End_year1992", "b_End_year1993","b_End_year1994","b_End_year1995","b_End_year1996") ### Need to change
  for(i in 1:7){#i=1 ## need to change
    pos.y<-(7:1)[i] ### need to change
    pos.x<-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"]
    lines(modoutput1[(modoutput1$term==rownameshere[i]),c("lower","upper")],rep(pos.y,2),col="darkgrey")
    points(pos.x,pos.y,cex=1.5,pch=19,col="darkblue")

    for(spsi in 1:length(spnum)){#
      pos.sps.i<-which(grepl(paste("[",spsi,"]",sep=""),mod.ranef$parameter,fixed=TRUE))
      jitt<-runif(1,0.05,0.5)
      pos.y.sps.i<-pos.y-jitt
      pos.x.sps.i<-mod.ranef[pos.sps.i[i],"mean"]
      
      lines(mod.ranef[pos.sps.i[i],c("25%","75%")],rep(pos.y.sps.i,2),
            col=alpha(my.pal[spsi], alphahere))
     
      points(pos.x.sps.i,pos.y.sps.i,cex=0.8,pch=19, col=alpha(my.pal[spsi], alphahere))
    }
    
  }
  par(xpd=TRUE) # so I can plot legend outside
  legend(leg1, leg2, sort(unique(gsub("_", " ", daters.fullrecords$species))), pch=rep(c(19,17),each=10), ### need to change
         col=alpha(my.pal[1:length(spnum)], alphahere),
         cex=1, bty="n", text.font=3)
  
}
