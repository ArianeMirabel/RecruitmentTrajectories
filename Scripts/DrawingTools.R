ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12);treatments<-list(T0,T1,T2,T3)

smoothTraj<-function(line){
  name<-names(line)
  line<-as.numeric(line)
  first<-(line[1]+line[2])/2;last<-(line[length(line)-1]+line[length(line)])/2
  ret<-c(first,unlist(lapply(2:(length(line)-1),function(step){return((line[step-1]+line[step]+line[step+1])/3)})),
         last)
  names(ret)<-name;return(ret)}

TrajectoryRec<-function(RecDB,s){
  invisible(lapply(RecDB,function(recind){
    Ylim=c(min(recind),max(recind)*s)
    plot(colnames(recind),recind[1,,"0.5"],type="n",ylim=Ylim,xlab="",ylab="")
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-recind[which(rownames(recind)%in%treatments[[tr]]),,]
      lapply(1:nrow(toplot),function(Li){
        lig<-t(toplot[Li,,]);lig<-lig[,which(!apply(lig,2,anyNA))]
        lig<-t(apply(lig,1,function(li){return(
          c(unlist(lapply(2:(length(li)-1),function(step){return((li[step-1]+li[step]+li[step+1])/3)})),li[length(li)]))}))
        lines(colnames(lig),lig["0.5",],col=ColorsTr[[tr]],lwd=2)
        polygon(c(colnames(lig),rev(colnames(lig))),c(lig["0.025",],rev(lig["0.975",])),
                col=rgb(0,0,0,alpha=0.05),border=NA)
      })}))
  }))}

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

PlotCWM<-function(TrajTraits){
  par(mfrow=c(2,dim(TrajTraits[[1]])[4]/2),mar=c(3,2,2,1),oma=c(2,2,1,1),no.readonly=TRUE)
for (Ntrait in 1:dim(TrajTraits[[1]])[4]){
  Toplot<-lapply(TrajTraits,function(tr){return(tr[,,,Ntrait])})
  Toplot<-lapply(Toplot,function(tr){return(tr[,which(colnames(tr)<=2015),])})
  Toplot<-lapply(Toplot,function(tr){
    ret<-do.call(rbind,lapply(c(0.025,0.5,0.975),function(quant){
      apply(apply(tr,c(1,2),function(rep){quantile(rep,probs=quant)}),2,mean)}))
    rownames(ret)<-c(0.025,0.5,0.975)
    colnames(ret)<-as.numeric(colnames(ret))-1984
    return(ret)})
  
  plot(colnames(Toplot[[1]]),Toplot[[1]]["0.5",],type="n",xaxt="n",xlab="",
       ylab="",cex.lab=1.5,ylim=c(min(unlist(Toplot)),max(unlist(Toplot))))
  axis(1,at=seq(5,30,5),labels=T) 
  mtext(c("Leaf Thickness","Chlorophyll Content", "Leaf Toughness","SLA","WD","Bark Thickness")[Ntrait],
        3,line = 0.5,cex=0.9)
  mtext(c("µm",expression(paste("g.",mm^-2,sep = "")), "N",expression(paste(mm^2,".",mg^-1,sep = "")),
          expression(paste("g.",cm^-3,sep = "")),"mm")[Ntrait],
        3,line = 0.5,adj=-0.05,cex=0.7)
  lapply(1:length(Toplot),function(tr){
    lines(colnames(Toplot[[tr]]),Toplot[[tr]]["0.5",],col=ColorsTr[tr],lwd=1.5)
    polygon(c(colnames(Toplot[[tr]]),rev(colnames(Toplot[[tr]]))),col=rgb(0,0,0,alpha=0.1),
            border=NA,c(Toplot[[tr]]["0.975",],rev(Toplot[[tr]]["0.025",])))
  })
}
mtext("Years since disturbance",side=1,at=0.90,cex=0.8,line=0.5,outer=TRUE)}

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

TrajectoryRec_fun<-function(RecDB_fun){
  Ylim=c(min(RecDB_fun,na.rm=T),max(RecDB_fun,na.rm=T))
  plot(colnames(RecDB_fun),RecDB_fun[1,,"0.5"],type="n",ylim=Ylim,xlab="",ylab="")
  invisible(lapply(1:length(treatments),function(tr){
    toplot<-RecDB_fun[which(rownames(RecDB_fun)%in%treatments[[tr]]),,]
    lapply(1:nrow(toplot),function(Li){
      lig<-t(toplot[Li,,]);lig<-lig[,which(!apply(lig,2,anyNA))]
      lig<-t(apply(lig,1,function(li){return(
        c(unlist(lapply(2:(length(li)-1),function(step){return((li[step-1]+li[step]+li[step+1])/3)})),li[length(li)]))}))
      lines(colnames(lig),lig["0.5",],col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(lig),rev(colnames(lig))),c(lig["0.025",],rev(lig["0.975",])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })}))
}

turnover<-function(TurnData){
  plot(colnames(TurnData),TurnData[1,,"0.5"],type="n",ylab="",xlab="",ylim=c(min(TurnData),max(TurnData)))
  invisible(lapply(1:length(treatments),function(tr){
    toplot<-TurnData[which(rownames(TurnData)%in%treatments[[tr]]),,]
    lapply(1:nrow(toplot),function(li){
      lines(colnames(toplot),smoothTraj(toplot[li,,"0.5"]),col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(toplot),rev(colnames(toplot))),
              c(smoothTraj(toplot[li,,"0.025"]),rev(smoothTraj(toplot[li,,"0.975"]))),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
  }))
}
