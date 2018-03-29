source("Scripts/Vernacular_handle.R")
load("DB/Alpha_Plots")
load("DB/Paracou_R_Subdivided_ok")
load("DB/BotanyGenus")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

#Plots of the different treatments
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3);names(treatments)<-c("0","1","2","3")

Recruits<-lapply(unique(Recruitment[,"campagne"]),function(y){
  return(Recruitment[which(Recruitment[,"campagne"]==y),])})
names(Recruits)<-unique(Recruitment[,"campagne"])
Recruits<-Recruits[order(names(Recruits))]
Recruits<-lapply(Recruits,function(yr){return(yr[which(yr[,"n_parcelle"]%in%1:12),])})
Recruits<-Recruits[which(lapply(Recruits,nrow)>0)]

dates<-names(Recruits[seq(2,29,2)])
Recruits3<-lapply(1:(length(dates)-1),function(y){
  return(do.call(rbind,Recruits[which(names(Recruits)>=dates[y] & names(Recruits)<dates[y+1])]))})
names(Recruits3)<-dates[-1]

#dates<-c("1985","1990","1995","2001","2006","2011","2017")

Nrep<-2
Nest<-lapply(1:Nrep,function(rep){
   Ret<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot),]
        ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        #ret<-as.character(yr[which(yr[,"n_parcelle"]==plot),"name"])
        #ret<-as.data.frame(Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
        #colnames(ret)<-"Species"
        #ret<-as.character(merge(ret, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
        ret<-tapply(ret,ret,length)
        
        reti<-LivingStand_all[["1984"]][[which(names(LivingStand_all[["1984"]])==plot)]]
        reti<-Replacement(reti,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
        #reti<-as.data.frame(Replacement(reti,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
        #colnames(reti)<-"Species"
        #reti<-as.character(merge(reti, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
        reti<-tapply(reti,reti,length)
        
        MetCom<-merge(ret,reti,by="row.names",all=T)
        rownames(MetCom)<-MetCom[,"Row.names"]; MetCom<-MetCom[,c("x","y")]
        MetCom[is.na(MetCom)]<-0
        
        nest<-merge(as.data.frame(apply(MetCom,1,min)),as.data.frame(apply(MetCom,1,max)),by="row.names")
        dimnames(nest)<-list(nest[,"Row.names"],c("rownames","min","max"));nest<-nest[,c("min","max")]
        nest<-merge(nest,MetCom,by="row.names")
        
        return((sum(nest[,"min"])+abs(sum(nest[,"x"])-sum(nest[,"y"])))/sum(nest[,"max"]))
    }))
  ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
  ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
  ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
  return(ret["ret"])
}))
return(Ret)})
Nest<-lapply(Nest,function(rep){return(smooth(rep,2))}) # moving average, path=2
Nest<-array(unlist(Nest),dim=c(nrow(Nest[[1]]),ncol(Nest[[1]]),Nrep),
                  dimnames=list(rownames(Nest[[1]]),as.numeric(names(Recruits3))-1984,1:Nrep))
Nest<-lapply(c(0.025,0.5,0.975),function(quant){return(
  apply(Nest,c(1,2),function(rep){return(quantile(rep,probs=quant))}))})
Nest<-array(unlist(Nest),dim=c(nrow(Nest[[1]]),ncol(Nest[[1]]),3),
                  dimnames=list(rownames(Nest[[1]]),colnames(Nest[[1]]),c(0.025,0.5,0.975)))

Turn<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
    ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
      ret<-yr[which(yr[,"n_parcelle"]==plot),]
      #ret<-as.character(yr[which(yr[,"n_parcelle"]==plot),"name"])
      ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
      #ret<-as.data.frame(Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
      #colnames(ret)<-"Species"
      #ret<-as.character(merge(ret, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
      ret<-tapply(ret,ret,length)
      
      reti<-LivingStand_all[["1991"]][[which(names(LivingStand_all[["1991"]])==plot)]]
      reti<-Replacement(reti,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
      #reti<-as.data.frame(Replacement(reti,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
      #colnames(reti)<-"Species"
      #reti<-as.character(merge(reti, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
      reti<-tapply(reti,reti,length)
      
      MetCom<-merge(ret,reti,by="row.names",all=T)
      rownames(MetCom)<-MetCom[,"Row.names"]; MetCom<-MetCom[,c("x","y")]
      MetCom[is.na(MetCom)]<-0
      
      return((sum(abs(MetCom[,"x"]-MetCom[,"y"]))-abs(sum(MetCom[,"x"])-sum(MetCom[,"y"])))/sum(apply(MetCom,1,max)))
    }))
    ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
    ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
    ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
    return(ret["ret"])
  }))
  return(Ret)})
Turn<-lapply(Turn,function(rep){return(smooth(rep,2))}) # moving average, path=2
Turn<-array(unlist(Turn),dim=c(nrow(Turn[[1]]),ncol(Turn[[1]]),Nrep),
                    dimnames=list(rownames(Turn[[1]]),as.numeric(names(Recruits3))-1984,1:Nrep))
Turn<-lapply(c(0.025,0.5,0.975),function(quant){return(
  apply(Turn,c(1,2),function(rep){return(quantile(rep,probs=quant))}))})
Turn<-array(unlist(Turn),dim=c(nrow(Turn[[1]]),ncol(Turn[[1]]),3),
                    dimnames=list(rownames(Turn[[1]]),colnames(Turn[[1]]),c(0.025,0.5,0.975)))


ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")
Toplot<-Turn
windows()
plot(colnames(Toplot),Toplot[1,,"0.5"],type="n",ylab="",xlab="",ylim=c(min(Toplot),max(Toplot)))
invisible(lapply(1:length(treatments),function(tr){
    toplot<-Toplot[which(rownames(Toplot)%in%treatments[[tr]]),,]
    lapply(1:nrow(toplot),function(li){
      lines(colnames(toplot),toplot[li,,"0.5"],col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(toplot),rev(colnames(toplot))),
              c(toplot[li,,"0.025"],rev(toplot[li,,"0.975"])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
}))
  
save(Nest,file="DB/Nestedness_toInit")
save(Turn,file="DB/Turnover_toInit")


  
