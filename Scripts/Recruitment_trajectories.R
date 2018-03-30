load("DB/Paracou_R_Subdivided_ok")
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

Nrep<-2
RecPun<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
  ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
  ret<-yr[which(yr[,"n_parcelle"]==plot),]
  ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
   return(expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id))
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
names(RecPun)<-c("Richness","Shannon","Simpson") 

RecAccum<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(1:length(Recruits3),function(Y){
      yr<-do.call(rbind,Recruits3[1:Y])
      if(length(which(duplicated(yr[,"i_arbre"])))!=0){print(Y)}
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot),]
        ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])#alphaGen)
        return(expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id))
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
names(RecAccum)<-c("Richness","Shannon","Simpson") 

save(RecAccum,file="DB/RecruitmentAccum")
save(RecPun,file="DB/RecruitmentPunctual")

########################################"
##### Null model

RecPun_DiffNull<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot),]
        ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        
        retN<-do.call(rbind,lapply(LivingStand_all[which(names(LivingStand_all)%in%unique(yr[,"campagne"]))],
                                   function(year){if(any(names(year)==plot)){return(year[[which(names(year)==plot)]])}}))
        retN<-retN[sample(rownames(retN),length(ret),replace=F),]
        retN<-Replacement(retN,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        retN<-expq(bcTsallis(as.AbdVector(tapply(retN,retN,length)),q=id,Correction="None"),q=id)
        ret<-expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id)
        return(ret-retN)
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
names(RecPun_DiffNull)<-c("Richness","Shannon","Simpson") 

RecAccum_DiffNull<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(1:length(Recruits3),function(Y){
      yr<-do.call(rbind,Recruits3[1:Y])
      if(length(which(duplicated(yr[,"i_arbre"])))!=0){print(Y)}
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot),]
        ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        
        retN<-do.call(rbind,lapply(LivingStand_all[which(names(LivingStand_all)%in%unique(yr[,"campagne"]))],
                                   function(year){if(any(names(year)==plot)){return(year[[which(names(year)==plot)]])}}))
        retN<-retN[sample(rownames(retN),length(ret),replace=F),]
        retN<-Replacement(retN,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        retN<-expq(bcTsallis(as.AbdVector(tapply(retN,retN,length)),q=id,Correction="None"),q=id)
        ret<-expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id)
        return(ret-retN)
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
names(RecAccum_DiffNull)<-c("Richness","Shannon","Simpson") 

save(RecPun_DiffNull,file="DB/RecruitmentPunctual_Nullmodel_Diff")
save(RecAccum_DiffNull,file="DB/RecruitmentAccum_Nullmodel_Diff")

RecPun_Null<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
       eff<-nrow(yr[which(yr[,"n_parcelle"]==plot),])
       ret<-do.call(rbind,lapply(LivingStand_all[which(names(LivingStand_all)%in%unique(yr[,"campagne"]))],
                    function(year){if(any(names(year)==plot)){return(year[[which(names(year)==plot)]])}}))
        ret<-ret[sample(rownames(ret),eff,replace=F),]
        ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        return(expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id))
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
names(RecPun_Null)<-c("Richness","Shannon","Simpson") 

RecAccum_Null<-lapply(0:2,function(id){
  recind<-lapply(1:Nrep,function(rep){
    Ret<-do.call(cbind,lapply(1:length(Recruits3),function(Y){
      yr<-do.call(rbind,Recruits3[1:Y])
      if(length(which(duplicated(yr[,"i_arbre"])))!=0){print(Y)}
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        eff<-nrow(yr[which(yr[,"n_parcelle"]==plot),])
        ret<-do.call(rbind,lapply(LivingStand_all[which(names(LivingStand_all)%in%unique(yr[,"campagne"]))],
                                  function(year){if(any(names(year)==plot)){return(year[[which(names(year)==plot)]])}}))
        ret<-ret[sample(rownames(ret),eff,replace=F),]
        ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        return(expq(bcTsallis(as.AbdVector(tapply(ret,ret,length)),q=id,Correction="None"),q=id))
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
names(RecAccum_Null)<-c("Richness","Shannon","Simpson") 

save(RecAccum_Null,file="RecruitmentAccum_Nullmodel")
save(RecPun_Null,file="RecruitmentPunctual_Nullmodel")

