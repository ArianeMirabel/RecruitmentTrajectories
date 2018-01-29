ColorsTr<-c("darkolivegreen","deepskyblue2","darkorange1","red2")

TrajectoryRec<-function(RecDB){
  Toplot<-lapply(c(1:4),function(t){
  toplot<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(RecDB[[t]][,,,q],c(1,2),function(x){return(quantile(x,probs = quant))}) )})
  toplot<-array(unlist(toplot),dim=c(nrow(toplot[[1]]),ncol(toplot[[1]]),3),
                dimnames=list(rownames(toplot[[1]]),
                              as.numeric(colnames(toplot[[1]]))-1984,c(0.025,0.5,0.975)))
  return(toplot)})

Ylim<-c(min(unlist(lapply(Toplot,function(tr){return(min(tr))}))),
        max(unlist(lapply(Toplot,function(tr){return(max(tr))}))))
plot(colnames(Toplot[[1]]),Toplot[[1]][1,,"0.5"],type="n",xaxt="n",xlab="",ylab="",
     ylim=Ylim)
axis(1,at=colnames(Toplot[[1]]),labels=T)

invisible(lapply(c(1:4),function(t){
  invisible(lapply(1:nrow(Toplot[[t]]),function(i){
    lines(colnames(Toplot[[t]]),Toplot[[t]][i,,"0.5"], col =  ColorsTr[t],lty = 1,lwd=1.5)
    polygon(c(colnames(Toplot[[t]]),rev(colnames(Toplot[[t]]))),
            c(Toplot[[t]][i,,"0.025"],rev(Toplot[[t]][i,,"0.975"])),
            col=rgb(0,0,0,alpha=0.1),border=NA)}))
}))}
