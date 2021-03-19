## Started 3 April 2019 ##
## By Cat - based off code from Lizzie's OSPREE plots ##

# Runs from models_stan_plotting.R #

# with COLORS for each species #

muplotfx2 <- function(modeloutput1,modeloutput2, mod.ran1,mod.ran2,nameforfig, width, height, ylim, xlim, leg1, leg2,leg3,leg4,name1,name2){
  spnum <- unique(dat$GEN.SPA) ### need to change
  #setEPS()
  #postscript(file.path(figpath, paste("", nameforfig, ".eps", sep="")),
   #          width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,5,3,0))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main="")
  axis(2, at=1:6, labels=rev(c("Chill", "Photo", "Force","Chill:\nPhoto","Chill:\nForce","Photo:\nForce")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  ppeffects <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  for(i in 1:6){#i=1 ## need to change
    pos.y<-(6:1)[i] ### need to change
    pos.x<-(modeloutput1[(modeloutput1$term==rownameshere[i]),"estimate"]-modeloutput2[(modeloutput2$term==rownameshere[i]),"estimate"])
    #pos.x1<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"])
    #pos.x2<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"])
    #lines(c(0,pos.x1),rep(pos.y,2),cex=10,lty=1,col="darkblue")
    lines(c(0,pos.x),rep(pos.y,2),lwd=6,lty=1,col="darkblue")
    #lines(c(0,pos.x2),rep(pos.y-.3,2),cex=10,lty=3,col="darkblue")
    #text(0,pos.y, labels="flowering",cex=.5,postion="right")
    #text(pos.x,pos.y, labels="budburst",cex=.5,position="left")
    #text(x=0,y=,pos.y-.15, labels="flowering",cex=.1,postion="right")
    #text(pos.x,pos.y-.15, labels="leafout",cex=.1,position="left")
    #text(0,pos.y-.3, labels="budburst",cex=.1,position="right")
    #text(pos.x,pos.y-.3, labels="leafout",cex=.1,position="left")
    values<-c()
    
    for(spsi in 1:length(spnum)){
      pos.sps.i<-which(grepl(paste("[",spsi,"]",sep=""),mod.ranef$parameter,fixed=TRUE))
      
      jitt<-0.075
      ii<-as.numeric(as.factor((spnum[spsi])))
      valueshere<-c(ii) 
      values<-sum(values,valueshere)
      pos.y.sps.i<-pos.y-jitt*values
      print(pos.y.sps.i)
      
      pos.x.sps.i<-mod.ran1[pos.sps.i[i],"mean"]-mod.ran2[pos.sps.i[i],"mean"]
      #points(pos.x.sps.i,pos.y.sps.i,cex=0.8,pch=17, col=alpha(my.pal[spsi], alphahere))
      lines(c(0,pos.x.sps.i),rep(pos.y.sps.i,2), lwd=2,col=alpha(my.pal[spsi], alphahere))

    }
    
  }
#  par(xpd=TRUE) # so I can plot legend outside
 # legend(leg1, leg2, sort(unique(gsub("_", " ", c("Comptonia perigrina","Viburnum acerifolium", "Ilex verticillata","Vaccinium corymbosum","Acer penyslyvanicum",
  #                                                "Prunus virginiana","Corylus cornuta","Prunus pensylvanica","Ilex mucronata","Acer rubrum")))), pch=rep(c(15),each=10), ### need to change
  #       col=alpha(my.pal[1:length(spnum)], alphahere),
   #      cex=.8, bty="n", text.font=3)
  #legend(leg3, leg4, legend=leg.txt, ### need to change
   #      cex=.8, bty="n", text.font=3)
  
  text(0,6.25, labels=name1)
  text(-20,6.25, labels=name2)
  text(20,6.25, labels=name2)
  }


muplotfx3 <- function(modeloutput1,modeloutput2, mod.ran1,mod.ran2,nameforfig, width, height, ylim, xlim, leg1, leg2,leg3,leg4,name1,name2){
  spnum <- unique(dat$GEN.SPA) ### need to change
  #setEPS()
  #postscript(file.path(figpath, paste("", nameforfig, ".eps", sep="")),
  #          width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,0,3,0))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main="")
  #axis(2, at=1:6) #labels=rev(c("Chill", "Photo", "Force","Chill:Photo","Chill:Force","Photo:Force")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  ppeffects <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  for(i in 1:6){#i=1 ## need to change
    pos.y<-(6:1)[i] ### need to change
    pos.x<-(modeloutput1[(modeloutput1$term==rownameshere[i]),"estimate"]-modeloutput2[(modeloutput2$term==rownameshere[i]),"estimate"])
    #pos.x1<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"])
    #pos.x2<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"])
    #lines(c(0,pos.x1),rep(pos.y,2),cex=10,lty=1,col="darkblue")
    lines(c(0,pos.x),rep(pos.y,2),lwd=6,lty=1,col="darkblue")
    #lines(c(0,pos.x2),rep(pos.y-.3,2),cex=10,lty=3,col="darkblue")
    #text(0,pos.y, labels="flowering",cex=.5,postion="right")
    #text(pos.x,pos.y, labels="budburst",cex=.5,position="left")
    #text(x=0,y=,pos.y-.15, labels="flowering",cex=.1,postion="right")
    #text(pos.x,pos.y-.15, labels="leafout",cex=.1,position="left")
    #text(0,pos.y-.3, labels="budburst",cex=.1,position="right")
    #text(pos.x,pos.y-.3, labels="leafout",cex=.1,position="left")
    values<-c()
    
    for(spsi in 1:length(spnum)){
      pos.sps.i<-which(grepl(paste("[",spsi,"]",sep=""),mod.ranef$parameter,fixed=TRUE))
      
      jitt<-0.075
      ii<-as.numeric(as.factor((spnum[spsi])))
      valueshere<-c(ii) 
      values<-sum(values,valueshere)
      pos.y.sps.i<-pos.y-jitt*values
      print(pos.y.sps.i)
      
      pos.x.sps.i<-mod.ran1[pos.sps.i[i],"mean"]-mod.ran2[pos.sps.i[i],"mean"]
      #points(pos.x.sps.i,pos.y.sps.i,cex=0.8,pch=17, col=alpha(my.pal[spsi], alphahere))
      lines(c(0,pos.x.sps.i),rep(pos.y.sps.i,2), lwd=2,col=alpha(my.pal[spsi], alphahere))
      
    }
    
  }
  #  par(xpd=TRUE) # so I can plot legend outside
  # legend(leg1, leg2, sort(unique(gsub("_", " ", c("Comptonia perigrina","Viburnum acerifolium", "Ilex verticillata","Vaccinium corymbosum","Acer penyslyvanicum",
  #                                                "Prunus virginiana","Corylus cornuta","Prunus pensylvanica","Ilex mucronata","Acer rubrum")))), pch=rep(c(15),each=10), ### need to change
  #       col=alpha(my.pal[1:length(spnum)], alphahere),
  #      cex=.8, bty="n", text.font=3)
  #legend(leg3, leg4, legend=leg.txt, ### need to change
  #      cex=.8, bty="n", text.font=3)
  
  text(0,6.25, labels=name1)
  text(-20,6.25, labels=name2)
  text(20,6.25, labels=name2)
}

muplotfx4 <- function(modeloutput1,modeloutput2, mod.ran1,mod.ran2,nameforfig, width, height, ylim, xlim, leg1, leg2,leg3,leg4,name1,name2){
  spnum <- unique(dat$GEN.SPA) ### need to change
  #setEPS()
  #postscript(file.path(figpath, paste("", nameforfig, ".eps", sep="")),
  #          width = width, height = height)
  par(xpd=FALSE)
  par(mar=c(5,0,3,8.5))
  plot(x=NULL,y=NULL, xlim=xlim, yaxt='n', ylim=ylim,
       xlab=xlab, ylab="", main="")
  #axis(2, at=1:6) #labels=rev(c("Chill", "Photo", "Force","Chill:Photo","Chill:Force","Photo:Force")), las=1) ### Need to change
  abline(v=0, lty=2, col="darkgrey")
  rownameshere <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  ppeffects <- c("b_Chill", "b_Light", "b_Force","b_Chill:Light","b_Chill:Force","b_Light:Force") ### Need to change
  for(i in 1:6){#i=1 ## need to change
    pos.y<-(6:1)[i] ### need to change
    pos.x<-(modeloutput1[(modeloutput1$term==rownameshere[i]),"estimate"]-modeloutput2[(modeloutput2$term==rownameshere[i]),"estimate"])
    #pos.x1<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput2[(modoutput2$term==rownameshere[i]),"estimate"])
    #pos.x2<-(modoutput3[(modoutput3$term==rownameshere[i]),"estimate"]-modoutput1[(modoutput1$term==rownameshere[i]),"estimate"])
    #lines(c(0,pos.x1),rep(pos.y,2),cex=10,lty=1,col="darkblue")
    lines(c(0,pos.x),rep(pos.y,2),lwd=6,lty=1,col="darkblue")
    #lines(c(0,pos.x2),rep(pos.y-.3,2),cex=10,lty=3,col="darkblue")
    #text(0,pos.y, labels="flowering",cex=.5,postion="right")
    #text(pos.x,pos.y, labels="budburst",cex=.5,position="left")
    #text(x=0,y=,pos.y-.15, labels="flowering",cex=.1,postion="right")
    #text(pos.x,pos.y-.15, labels="leafout",cex=.1,position="left")
    #text(0,pos.y-.3, labels="budburst",cex=.1,position="right")
    #text(pos.x,pos.y-.3, labels="leafout",cex=.1,position="left")
    values<-c()
    
    for(spsi in 1:length(spnum)){
      pos.sps.i<-which(grepl(paste("[",spsi,"]",sep=""),mod.ranef$parameter,fixed=TRUE))
      
      jitt<-0.075
      ii<-as.numeric(as.factor((spnum[spsi])))
      valueshere<-c(ii) 
      values<-sum(values,valueshere)
      pos.y.sps.i<-pos.y-jitt*values
      print(pos.y.sps.i)
      
      pos.x.sps.i<-mod.ran1[pos.sps.i[i],"mean"]-mod.ran2[pos.sps.i[i],"mean"]
      #points(pos.x.sps.i,pos.y.sps.i,cex=0.8,pch=17, col=alpha(my.pal[spsi], alphahere))
      lines(c(0,pos.x.sps.i),rep(pos.y.sps.i,2), lwd=2,col=alpha(my.pal[spsi], alphahere))
      
    }
    
  }
  
  text(0,6.25, labels=name1)
  text(-20,6.25, labels=name2)
  text(20,6.25, labels=name2)
  
  par(xpd=TRUE) # so I can plot legend outside
  legend(leg1, leg2, sort(unique(gsub("_", " ", c("Comptonia perigrina","Viburnum acerifolium", "Ilex verticillata","Vaccinium corymbosum","Acer penyslyvanicum",
                                                  "Prunus virginiana","Corylus cornuta","Prunus pensylvanica","Ilex mucronata","Acer rubrum")))), pch=rep(c(15),each=10), ### need to change
         col=alpha(my.pal[1:length(spnum)], alphahere),
         cex=.8, bty="n", text.font=3)  
}
