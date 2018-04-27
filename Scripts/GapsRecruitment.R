Rec<-read.table("DB/GapsRecruitment.csv",sep=";",header=T)

load("DB/Paracou_R_Subdivided_ok")

Recruitment<-merge(Recruitment,Rec[,c("i_arbre","PresenceExploit")],by="i_arbre",all=T)
Recruitment[which(is.na(Recruitment[,"PresenceExploit"])),"PresenceExploit"]<-0
colnames(Recruitment)<-c("i_arbre","n_parcelle","nomPilote","code_vivant","campagne","Famille","Genre","Espece","name","gap")

load("DB/Alpha_Plots")
source("Scripts/Vernacular_handle.R")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

#Plots of the different treatments
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

Recruits<-lapply(unique(Recruitment[,"campagne"]),function(y){
  return(Recruitment[which(Recruitment[,"campagne"]==y),])})
names(Recruits)<-unique(Recruitment[,"campagne"])
Recruits<-Recruits[order(names(Recruits))]
Recruits<-lapply(Recruits,function(yr){return(yr[which(yr[,"n_parcelle"]%in%1:15),])})
Recruits<-Recruits[which(lapply(Recruits,nrow)>0)]

dates<-names(Recruits[seq(2,29,2)])
Recruits3<-lapply(1:(length(dates)-1),function(y){
  return(do.call(rbind,Recruits[which(names(Recruits)>=dates[y] & names(Recruits)<dates[y+1])]))})
names(Recruits3)<-dates[-1]

Nrep<-50
RecG<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot & yr[,"gap"]==1),]
        if(nrow(ret)!=0){ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        return(expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id))}
        else(return(0))
      }))
      ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
      ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
      ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
      return(ret["ret"])
    }))
    return(Ret)})
  recind<-lapply(recind,function(rep){return(smooth(rep,2))}) # moving average, path=2
  recind<-array(unlist(recind),dim=c(nrow(recind[[1]]),ncol(recind[[1]]),Nrep),
                dimnames=list(rownames(recind[[1]]),names(Recruits3),1:Nrep))
  recind<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(recind,c(1,2),function(rep){return(quantile(rep,probs = quant,na.rm=T))}))
  })
  recind<-array(unlist(recind),dim=c(nrow(recind[[1]]),ncol(recind[[1]]),3),
                dimnames=list(rownames(recind[[1]]),as.numeric(colnames(recind[[1]]))-1984,c(0.025,0.5,0.975)))
  return(recind)
})
names(RecG)<-c("Richness","Shannon","Simpson") 

RecNG<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot & yr[,"gap"]==0),]
        if(nrow(ret)!=0){ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        return(expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id))}
        else(return(0))
      }))
      ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
      ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
      ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
      return(ret["ret"])
    }))
    return(Ret)})
  recind<-lapply(recind,function(rep){return(smooth(rep,2))}) # moving average, path=2
  recind<-array(unlist(recind),dim=c(nrow(recind[[1]]),ncol(recind[[1]]),Nrep),
                dimnames=list(rownames(recind[[1]]),names(Recruits3),1:Nrep))
  recind<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(recind,c(1,2),function(rep){return(quantile(rep,probs = quant,na.rm=T))}))
  })
  recind<-array(unlist(recind),dim=c(nrow(recind[[1]]),ncol(recind[[1]]),3),
                dimnames=list(rownames(recind[[1]]),as.numeric(colnames(recind[[1]]))-1984,c(0.025,0.5,0.975)))
  return(recind)
})
names(RecNG)<-c("Richness","Shannon","Simpson") 

save(RecG,file="DB/RecruitmentGapCorr")
save(RecNG,file="DB/RecruitmentNonGapCorr")

gaps<-function(RecDBg,RecDBng){
  invisible(lapply(1:length(RecDB),function(ind){
    recindg<-RecDBg[[ind]]
    recindng<-RecDBng[[ind]]
    Ylim=c(min(recindg,recindng),max(recindg,recindng))
    plot(colnames(recindg),recindg[1,,"0.5"],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n",main=names(RecDB)[ind])
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-recindg[which(rownames(recindg)%in%treatments[[tr]]),,]
      lapply(1:nrow(toplot),function(Li){
        lines(colnames(toplot),toplot[Li,,"0.5"],col=ColorsTr[[tr]],lwd=1.5)
        polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[Li,,"0.025"],rev(toplot[Li,,"0.975"])),
                col=rgb(0,0,0,alpha=0.05),border=NA)
      })
    }))
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-recindng[which(rownames(recindng)%in%treatments[[tr]]),,]
      lapply(1:nrow(toplot),function(Li){
        lines(colnames(toplot),toplot[Li,,"0.5"],col=ColorsTr[[tr]],lwd=1.5,lty=3)
        polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[Li,,"0.025"],rev(toplot[Li,,"0.975"])),
                col=rgb(0,0,0,alpha=0.05),border=NA)})}))
}))

}

windows(width=60,height=20)
par(mfcol=c(1,3),oma=c(1,1,1,0),mar=c(2,2,2,0),no.readonly=TRUE)
gaps(RecNG,RecNG)


