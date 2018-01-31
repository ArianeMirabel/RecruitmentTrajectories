ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12);treatments<-list(T0,T1,T2,T3)

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

TrajectoryRec_vsNull<-function(RecDB_obs,RecDb_null){
  par(mfrow=c(4,3),mar=c(2,3,1,1),oma=c(2,2,2.5,1),no.readonly=TRUE)
  layout(matrix(1:12, 4, 3, byrow = FALSE))
for (ind in c("Richness","Shannon","Simpson")){
  Toplot_Null<-lapply(RecDb_null,function(tr){return(tr[,,,ind])})
  Toplot_Null<-lapply(RecDb_null,function(tr){
    colnames(tr)<-as.numeric(colnames(tr))-1984
    ret<-lapply(c(0.025,0.5,0.975),function(quant){
      apply(tr,c(1,2),function(x){return(quantile(x,probs=quant))})})
    names(ret)<-c(0.025,0.5,0.975)
    return(ret)})
  
  Toplot<-lapply(RecDB_obs,function(tr){return(tr[,,,ind])})
  Toplot<-lapply(Toplot,function(tr){
    colnames(tr)<-as.numeric(colnames(tr))-1984
    ret<-lapply(c(0.025,0.5,0.975),function(quant){
      apply(tr,c(1,2),function(x){return(quantile(x,probs=quant))})})
    names(ret)<-c(0.025,0.5,0.975)
    return(ret)})
  
  lapply(1:length(Toplot_Null),function(tr){
    plot(colnames(Toplot_Null[[tr]][[1]]),Toplot_Null[[tr]][["0.5"]][1,],type="n",xaxt="n",xlab="",
         ylab="",cex.lab=1.5,
         ylim=c(min(unlist(Toplot_Null[[tr]]),unlist(Toplot[[tr]]),na.rm=T),
                max(unlist(Toplot_Null[[tr]]),unlist(Toplot[[tr]]),na.rm=T)))
    axis(1,at=seq(5,30,5),labels=T) 
    toplot_Null<-Toplot_Null[[tr]][["0.5"]]
    toplot<-Toplot[[tr]][["0.5"]]
    invisible(lapply(1:nrow(toplot_Null),function(li){
      lines(colnames(toplot_Null),toplot_Null[li,],col=ColorsTr[tr],lwd=1.5,lty=2)
      #polygon(c(colnames(toplot_Null),rev(colnames(toplot_Null))),col=rgb(0,0,0,alpha=0.1),border=NA,
      #c(Toplot_Null[[tr]][["0.975"]][li,],rev(Toplot_Null[[tr]][["0.025"]][li,])))
      
      lines(colnames(toplot),toplot[li,],col=ColorsTr[tr],lwd=2)
      polygon(c(colnames(toplot),rev(colnames(toplot))),col=rgb(0,0,0,alpha=0.1),border=NA,
              c(Toplot[[tr]][["0.975"]][li,],rev(Toplot[[tr]][["0.025"]][li,])))
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

turnover<-function(TurnData){
  plot(colnames(TurnData[[1]]),TurnData[[1]][1,],type="n",xaxt="n",xlab="",ylab="",
       ylim=c(min(unlist(lapply(TurnData,min))),max(unlist(lapply(TurnData,max)))))
  axis(1,at=seq(5,30,5),labels=TRUE)    
  
  invisible(lapply(1:4,function(t){  
    
    toplot<-TurnData[[t]]
    
    lines(colnames(toplot),toplot["0.5",], col =  ColorsTr[t],lty = 1,lwd=3)
    polygon(c(colnames(toplot),rev(colnames(toplot))),
            c(toplot["0.975",],rev(toplot["0.025",])),
            col=rgb(0,0,0,alpha=0.1),border=NA)
  }))
  mtext("Recruitment turnover compared to current stand",adj=0,line=1,cex=1.2)
  mtext("Years since disturbance",side=1,line=2,adj=1)
}
