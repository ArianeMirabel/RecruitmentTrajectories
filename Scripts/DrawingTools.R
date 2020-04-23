load("DB/LostAGB")
norm<-(AGBloss_cor[,"AGB"]-min(AGBloss_cor[,"AGB"]))/(max(AGBloss_cor[,"AGB"])-min(AGBloss_cor[,"AGB"]))
ColorsDist <- colorRampPalette(c("darkolivegreen2","gold","orangered","darkred"))(12)[as.numeric(cut(norm, breaks = 12))]


TrajectoryDiffNull<-function(RecDB,RecDB_Diff){
  invisible(lapply(c(1,3),function(ind){
    recind<-RecDB[[ind]]
    Ylim=c(min(recind),max(recind))
    par(mar=c(0, 2, 3, 1))
    plot(colnames(recind),recind[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n")
    mtext(c("(a) Taxonomic Richness", "", "(b) Taxonomic Evenness")[ind],side=3,line=3.5,cex=0.9)
    mtext(c("Equivalent\ndiversity", "", "")[ind],side=3,line=0.5,cex=0.77,adj=0)
    lapply(1:nrow(recind),function(Li){
      lines(colnames(recind),recind[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2)
      polygon(c(colnames(recind),rev(colnames(recind))),
              c(recind[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(recind[AGBloss_cor[,"plot"],,"0.975"][Li,])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
    par(mar=c(1, 2, 0, 1))
    recindN<-RecDB_Diff[[ind]]
    Ylim=c(min(recindN),max(recindN))
    plot(colnames(recindN),recind[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n")
    axis(1,at=colnames(recindN),labels=T) 
    abline(h=0,lty=1,lwd=1)
    lapply(1:nrow(recindN),function(Li){
      lines(colnames(recindN),recindN[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2,lty=2)
      polygon(c(colnames(recindN),rev(colnames(recindN))),
              c(recindN[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(recindN[AGBloss_cor[,"plot"],,"0.975"][Li,])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
  }))
  mtext(expression(paste("Communities diversity, ",'H'['obs']^q)),side=2,line=0.5,adj=0,at=0.45,cex=0.8,outer=TRUE)
  mtext(expression(paste('H'['obs']^q," - ",'H'['null']^q)),
        side=2,line=0.5,adj=0,at=0.1,cex=0.9,outer=TRUE)
}

TrajectoryRec_fun<-function(RecDB_fun,RecDB_fun_Diff){
  Ylim=c(min(RecDB_fun,na.rm=T),max(RecDB_fun,na.rm=T))
  par(mar=c(0, 2, 3, 1))
  plot(colnames(RecDB_fun),RecDB_fun[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n")
  mtext("(c) Functional Diversity",side=3,line=3.5,cex=0.9)
  lapply(1:nrow(RecDB_fun),function(Li){
    lines(colnames(RecDB_fun),RecDB_fun[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2)
    polygon(c(colnames(RecDB_fun),rev(colnames(RecDB_fun))),
            c(RecDB_fun[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(RecDB_fun[AGBloss_cor[,"plot"],,"0.975"][Li,])),
            col=rgb(0,0,0,alpha=0.05),border=NA)
  })
  par(mar=c(1, 2, 0, 1))
  Ylim=c(min(RecDB_fun_Diff,na.rm=T),max(RecDB_fun_Diff,na.rm=T))
  plot(colnames(RecDB_fun_Diff),RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="")
  abline(h=0,lty=1,lwd=1)
  lapply(1:nrow(RecDB_fun_Diff),function(Li){
    lines(colnames(RecDB_fun_Diff),RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2,lty=2)
    polygon(c(colnames(RecDB_fun_Diff),rev(colnames(RecDB_fun_Diff))),
            c(RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.975"][Li,])),
            col=rgb(0,0,0,alpha=0.03),border=NA)
  })
  mtext("Years since disturbance",side=1,at=0.85,line=1.5,outer=TRUE)
}


PlotCWM<-function(TrajTraits){
  for (Ntrait in 1:length(TrajTraits)){
    toplot<-Rec_CWM[[Ntrait]][,,"0.5"]
    nm<-colnames(toplot)
    toplot<-as.data.frame(cbind(rownames(toplot),toplot))
    
    toplot<-reshape(toplot, direction="long", varying = 2:14, idvar="V1", v.names="val", timevar="time", times=nm)
    toplot<-as.data.frame(apply(toplot,2,as.numeric))
    plot(x=toplot[,"time"],y=toplot[,"val"],type="n",xlab="time",ylab=names(Rec_CWM)[Ntrait])
    
    lapply(1:nrow(AGBloss_cor),function(Li){
      dat<-toplot[which(toplot[,"V1"]==AGBloss_cor[Li,"plot"]),]
      points(x=dat[,"time"],y=dat[,"val"],col=ColorsDist[Li],pch=20,cex=0.5)
      #gam_mod <- gam(val ~ s(time,k=5), data = dat, sp=0.1)
      gam_mod <- gam(val ~ s(time,k=5), data = dat, method="REML")
      xseq<-seq(from=5,to=30,length.out = 100)
      pred<-predict.gam(gam_mod,data.frame(time=xseq),type="response",se.fit=T)
      lines(y=pred$fit, x=xseq,col=ColorsDist[Li],lwd=1.5)
      polygon(c(xseq,rev(xseq)),c(pred$se.fit,rev(pred$se.fit)),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
  }
}

legendCWM<-function(){
  mtext("Leaf thickness\n",at=0.13,line=-1.5,outer=TRUE,cex=0.9)
  mtext(expression(paste(mu, "m",sep = "")),at=0.08,line=-1.9,outer=TRUE,cex=0.9)
  mtext("Leaf chlorophyll content\n",at=0.4,line=-1.5,outer=TRUE,cex=0.9)
  mtext(expression(paste("g.",mm^-2,sep = "")),at=0.34,line=-1.9,outer=TRUE,cex=0.9)
  mtext("Leaf toughness\n",at=0.64,line=-1.5,outer=TRUE,cex=0.9)
  mtext("N",at=0.56,line=-1.9,outer=TRUE,cex=0.9)
  mtext("SLA\n",at=0.88,line=-1.5,outer=TRUE,cex=0.9)
  mtext(expression(paste(mm^2,".",mg^-1,sep = "")),at=0.84,line=-1.9,outer=TRUE,cex=0.9)
  
  mtext("WSG\n",at=0.13,line=-16.3,outer=TRUE,cex=0.9)
  mtext(expression(paste("g.",cm^-3,sep = "")),at=0.08,line=-16.5,outer=TRUE,cex=0.9)
  mtext("Bark thickness\n",at=0.4,line=-16.6,outer=TRUE,cex=0.9)
  mtext("mm",at=0.32,line=-16.2,outer=TRUE,cex=0.9)
  mtext("Hmax\n",at=0.64,line=-16.2,outer=TRUE,cex=0.9)
  mtext("m",at=0.56,line=-16.2,outer=TRUE,cex=0.9)
}

FDiversity<-function(FdivDB){
  plot(colnames(FdivDB[[1]]),FdivDB[[1]][1,],type="n",xaxt="n",xlab="",
       ylab="",cex.lab=1.5,ylim=c(min(unlist(FdivDB)),max(unlist(FdivDB))))
  axis(1,at=seq(5,30,5),labels=T) 
  invisible(lapply(1:length(FdivDB),function(tr){
    lines(colnames(FdivDB[[tr]]),FdivDB[[tr]]["0.5",],col=ColorsTr[tr],lwd=2)
    
    polygon(c(colnames(FdivDB[[tr]]),rev(colnames(FdivDB[[tr]]))),
            col=rgb(0,0,0,alpha=0.1),border=NA,
            c(FdivDB[[tr]]["0.975",],rev(FdivDB[[tr]]["0.025",])))
  }))}

TrajectoryRec_fun<-function(RecDB_fun,RecDB_fun_Diff){
  Ylim=c(min(RecDB_fun,na.rm=T),max(RecDB_fun,na.rm=T))
  par(mar=c(0, 2, 3, 1))
  plot(colnames(RecDB_fun),RecDB_fun[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n")
  mtext("(c) Functional Diversity",side=3,line=3.5,cex=0.9)
  lapply(1:nrow(RecDB_fun),function(Li){
    lines(colnames(RecDB_fun),RecDB_fun[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2)
    polygon(c(colnames(RecDB_fun),rev(colnames(RecDB_fun))),
            c(RecDB_fun[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(RecDB_fun[AGBloss_cor[,"plot"],,"0.975"][Li,])),
            col=rgb(0,0,0,alpha=0.05),border=NA)
  })
  par(mar=c(1, 2, 0, 1))
  Ylim=c(min(RecDB_fun_Diff,na.rm=T),max(RecDB_fun_Diff,na.rm=T))
  plot(colnames(RecDB_fun_Diff),RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="")
  abline(h=0,lty=1,lwd=1)
  lapply(1:nrow(RecDB_fun_Diff),function(Li){
    lines(colnames(RecDB_fun_Diff),RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2,lty=2)
    polygon(c(colnames(RecDB_fun_Diff),rev(colnames(RecDB_fun_Diff))),
            c(RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(RecDB_fun_Diff[AGBloss_cor[,"plot"],,"0.975"][Li,])),
            col=rgb(0,0,0,alpha=0.03),border=NA)
  })
  mtext("Years since disturbance",side=1,at=0.85,line=1.5,outer=TRUE)
}

turnover<-function(TurnData){
  plot(colnames(TurnData),TurnData[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylab="",xlab="",ylim=c(min(TurnData),max(TurnData)),bty="n")
  invisible(lapply(1:nrow(TurnData),function(li){
      lines(colnames(TurnData),TurnData[AGBloss_cor[,"plot"],,"0.5"][li,],col=ColorsDist[li],lwd=2)
      polygon(c(colnames(TurnData),rev(colnames(TurnData))),
              c(TurnData[AGBloss_cor[,"plot"],,"0.025"][li,],rev(TurnData[AGBloss_cor[,"plot"],,"0.975"][li,])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    }))
}

Gaps<-function(RecDBg,RecDBng){
  par(mfrow=c(2,3),oma=c(3,4,4,0),no.readonly=TRUE)
  par(mar=c(0, 2, 2, 1))
  invisible(lapply(1:length(RecDBng),function(ind){
    recind<-RecDBng[[ind]]
    Ylim=c(min(recind),max(recind))
    plot(colnames(recind),recind[1,,"0.5"],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n",main=names(RecDBng)[ind])
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-recind[which(rownames(recind)%in%treatments[[tr]]),,]
      lapply(1:nrow(toplot),function(Li){
        lines(colnames(toplot),toplot[Li,,"0.5"],col=ColorsTr[[tr]],lwd=1.5)
        polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[Li,,"0.025"],rev(toplot[Li,,"0.975"])),
                col=rgb(0,0,0,alpha=0.05),border=NA)
      })
    }))
  }))
  
  par(mar=c(1, 2, 0, 1))
  invisible(lapply(1:length(RecDBg),function(ind){
    recind<-RecDBg[[ind]]
    Ylim=c(min(recind),max(recind))
    plot(colnames(recind),recind[1,,"0.5"],type="n",ylim=Ylim,xlab="",ylab="")
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-recind[which(rownames(recind)%in%treatments[[tr]]),,]
      lapply(1:nrow(toplot),function(Li){
        lines(colnames(toplot),toplot[Li,,"0.5"],col=ColorsTr[[tr]],lwd=1.5)
        polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[Li,,"0.025"],rev(toplot[Li,,"0.975"])),
                col=rgb(0,0,0,alpha=0.05),border=NA)
      })
    }))
  }))
  mtext("Equivalent\ndiversity",side=3,adj=0,line=-1,cex=1.2,outer=TRUE)
  mtext("Years since disturbance",side=1,at=0.85,line=2,cex=1.2,outer=TRUE)
  mtext("Undisturbed",side=2,line=1,cex=1.2,at=0.75,outer=TRUE)
  mtext("Gaps",side=2,line=1,cex=1.2,at=0.25,outer=TRUE)
}


##############################################################"


TrajectoryRec_vsNull<-function(RecDB_obs,RecDB_null){
  par(mfrow=c(4,3),mar=c(2,3,1,1),oma=c(2,2,2.5,1),no.readonly=TRUE)
  layout(matrix(1:12, 4, 3, byrow = FALSE))
  for (ind in c("Richness","Shannon","Simpson")){
    lapply(1:length(treatments),function(tr){
      toplot<-RecDB_obs[[ind]][which(rownames(RecDB_obs[[ind]])%in%treatments[[tr]]),,]
      toplotN<-RecDB_null[[ind]][which(rownames(RecDB_null[[ind]])%in%treatments[[tr]]),,]
      plot(colnames(toplotN),toplotN[1,,"0.5"],type="n",xaxt="n",xlab="",
           ylab="",cex.lab=1.5,ylim=c(min(min(toplot),min(toplotN)),max(max(toplot),max(toplotN))))
      axis(1,at=seq(5,30,5),labels=T) 
      invisible(lapply(1:nrow(toplot),function(li){
        lines(colnames(toplotN),toplotN[li,,"0.5"],col=ColorsTr[tr],lwd=1.5,lty=2)
        polygon(c(colnames(toplotN),rev(colnames(toplotN))),col=rgb(0,0,0,alpha=0.05),border=NA,
                c(toplotN[li,,"0.025"],rev(toplotN[li,,"0.975"])))
        
        lines(colnames(toplot),smoothTraj(toplot[li,,"0.5"]),col=ColorsTr[tr],lwd=2)
        polygon(c(colnames(toplot),rev(colnames(toplot))),col=rgb(0,0,0,alpha=0.1),border=NA,
                c(smoothTraj(toplot[li,,"0.025"]),rev(smoothTraj(toplot[li,,"0.975"]))))
      }))
    })
  }}


PlotCWM_old<-function(TrajTraits){
  for (Ntrait in 1:length(TrajTraits)){
    Toplot<-TrajTraits[[Ntrait]]
    plot(colnames(Toplot),Toplot[1,,"0.5"],type="n",xaxt="n",xlab="",
         ylab="",cex.lab=1.5,ylim=c(min(Toplot),max(Toplot)))
    axis(1,at=seq(5,30,5),labels=T) 
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-Toplot[which(rownames(Toplot)%in%treatments[[tr]]),,]  
      lapply(1:nrow(toplot),function(Li){
        lines(colnames(toplot),toplot[Li,,"0.5"],col=ColorsTr[[tr]],lwd=2)
        polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[Li,,"0.025"],rev(toplot[Li,,"0.975"])),
                col=rgb(0,0,0,alpha=0.05),border=NA)})
    }))
  }
}


### Main figure

.f = function() {
  
  load("DB/Turnover_toInit")
  load("DB/RecruitmentPunctual");load("DB/RecruitmentPunctual_Nullmodel_Diff")
  load("DB/RecruitmentPunctual_Functional");load("DB/RecruitmentPunctual_Functional_Nullmodel_Diff")
  load("DB/LostAGB")
  
  RecDB <-RecPun
  RecDB_fun <-RecPun_Fun
  
  #windows(width=100, height=40)
  par(mfcol=c(1,4),oma=c(5,2,5,3),no.readonly=TRUE)
  invisible(lapply(c(1,3),function(ind){
    recind<-RecDB[[ind]]
    Ylim=c(min(recind),max(recind))
    par(mar=c(0, 2, 3, 1))
    plot(colnames(recind),recind[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="")
    mtext(c("(a) Taxonomic Richness", "", "(b) Taxonomic Evenness")[ind],side=3,line=3.5,cex=0.9)
    mtext("Equivalent\ndiversity",side=3,line=0.5,cex=0.77,adj=0)
    lapply(1:nrow(recind),function(Li){
      lines(colnames(recind),recind[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2)
      polygon(c(colnames(recind),rev(colnames(recind))),
              c(recind[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(recind[AGBloss_cor[,"plot"],,"0.975"][Li,])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
  }))
  
  Ylim=c(min(RecDB_fun,na.rm=T),max(RecDB_fun,na.rm=T))
  par(mar=c(0, 2, 3, 1))
  plot(colnames(RecDB_fun),RecDB_fun[AGBloss_cor[,"plot"],,"0.5"][1,],type="n",ylim=Ylim,xlab="",ylab="")
  mtext("(c) Functional Diversity",side=3,line=3.5,cex=0.9)
  mtext("Equivalent\ndiversity",side=3,line=0.5,cex=0.77,adj=0)
  lapply(1:nrow(RecDB_fun),function(Li){
    lines(colnames(RecDB_fun),RecDB_fun[AGBloss_cor[,"plot"],,"0.5"][Li,],col=ColorsDist[Li],lwd=2)
    polygon(c(colnames(RecDB_fun),rev(colnames(RecDB_fun))),
            c(RecDB_fun[AGBloss_cor[,"plot"],,"0.025"][Li,],rev(RecDB_fun[AGBloss_cor[,"plot"],,"0.975"][Li,])),
            col=rgb(0,0,0,alpha=0.05),border=NA)
  })
  
  turnover(Turn)  
  mtext(c("(d) Taxonomic turnover"),side=3,line=3.5,cex=0.9)
  mtext("Years since disturbance",side=1,adj=1,outer=T,line=4,cex=1.1)
  mtext("Post-disturbance Recruitment Trajectories",side=3,outer=T,line=3.5,cex=1.1)
  
  legd<-rep(NA, 12);legd[c(1,6,12)]<-round(sort(AGBloss_cor[,"AGB"])[c(1,6,12)])
  legend("right",inset=c(-0.3,0),xpd=NA,legend=rev(legd),fill=rev(ColorsDist),bty="n",title="Lost AGB\n(%)\n", y.intersp = 0.5, border=NA)
  legend(x=31,y=0.05,xpd=NA,legend = "CI 95%",fill=rgb(0,0,0,alpha=0.1),border=NA,bty="n",x.intersp = 0.4)
  
}


















