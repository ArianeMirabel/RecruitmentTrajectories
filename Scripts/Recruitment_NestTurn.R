library(RColorBrewer)
setwd("P:/Private/Taff/These/Codes/DemographicDynamics")

source("Vernacular_handle.R")
source("Overview_functions.R")

load("DB/Alpha_Plots")
load("DB/Paracou_R_Subdivided_ok")
load("DB/BotanyGenus")

#Plots of the different treatments
T0<-c(1,6,11) 
######################  Warning! Warning! Pour l'instant je n'ai pas mis la 16, par ce que ça fout le bordel
T1<-c(2,7,9)
T2<-c(3,5,10)
T3<-c(4,8,12)

treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("0","1","2","3")

Recruits<-lapply(unique(Recruitment[,"campagne"]),function(y){
  return(Recruitment[which(Recruitment[,"campagne"]==y),])})
names(Recruits)<-unique(Recruitment[,"campagne"])
Recruits<-Recruits[order(names(Recruits))]
Recruits<-lapply(Recruits,function(yr){return(yr[which(yr[,"n_parcelle"]%in%1:15),])})
Recruits<-Recruits[which(lapply(Recruits,nrow)>0)]

dates<-names(Recruits)
sequence<-seq(1,29,3)
Recruits3<-lapply(1:(length(sequence)-1),function(i){return(do.call(rbind,
                                                                    Recruits[which(names(Recruits)>=dates[sequence[i]] & names(Recruits)<dates[sequence[i+1]])]))})
names(Recruits3)<-dates[sequence[-1]]
Recruits3<-lapply(Recruits3,function(yr){return(yr[which(yr[,"n_parcelle"]%in%1:12),])})

ColorsTr<-c("chartreuse3","deepskyblue2","darkorange1","red2")

dates<-c("1985","1990","1995","2001","2006","2011","2017")

Nrep<-50
Nest<-lapply(1:Nrep,function(rep){
   Ret<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        ret<-yr[which(yr[,"n_parcelle"]==plot),]
        #ret<-as.character(yr[which(yr[,"n_parcelle"]==plot),"name"])
        ret<-as.data.frame(Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
        colnames(ret)<-"Species"
        ret<-as.character(merge(ret, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
        ret<-tapply(ret,ret,length)
        
        reti<-LivingStand_all[["1984"]][[which(names(LivingStand_all[["1984"]])==plot)]]
        reti<-as.data.frame(Replacement(reti,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
        colnames(reti)<-"Species"
        reti<-as.character(merge(reti, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
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
colnames(Ret)<-as.numeric(names(Recruits3))-1984
return(Ret)})
Nest<-array(unlist(Nest),dim=c(nrow(Nest[[1]]),ncol(Nest[[1]]),Nrep),
                  dimnames=list(rownames(Nest[[1]]),colnames(Nest[[1]]),1:Nrep))
Nest<-lapply(c(0.025,0.5,0.975),function(quant){return(
  apply(Nest,c(1,2),function(rep){return(quantile(rep,probs=quant))}))})
Nest<-array(unlist(Nest),dim=c(nrow(Nest[[1]]),ncol(Nest[[1]]),3),
                  dimnames=list(rownames(Nest[[1]]),colnames(Nest[[1]]),c(0.025,0.5,0.975)))

Turn<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
    ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
      ret<-yr[which(yr[,"n_parcelle"]==plot),]
      #ret<-as.character(yr[which(yr[,"n_parcelle"]==plot),"name"])
      ret<-as.data.frame(Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
      colnames(ret)<-"Species"
      ret<-as.character(merge(ret, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
      ret<-tapply(ret,ret,length)
      
      reti<-LivingStand_all[["1991"]][[which(names(LivingStand_all[["1991"]])==plot)]]
      reti<-as.data.frame(Replacement(reti,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]]))
      colnames(reti)<-"Species"
      reti<-as.character(merge(reti, RefBota,by.x="Species",by.y="row.names",all.x=T)[,"Genre"])
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
  colnames(Ret)<-as.numeric(names(Recruits3))-1984
  return(Ret)})
Turn<-array(unlist(Turn),dim=c(nrow(Turn[[1]]),ncol(Turn[[1]]),Nrep),
                    dimnames=list(rownames(Turn[[1]]),colnames(Turn[[1]]),1:Nrep))
Turn<-lapply(c(0.025,0.5,0.975),function(quant){return(
  apply(Turn,c(1,2),function(rep){return(quantile(rep,probs=quant))}))})
Turn<-array(unlist(Turn),dim=c(nrow(Turn[[1]]),ncol(Turn[[1]]),3),
                    dimnames=list(rownames(Turn[[1]]),colnames(Turn[[1]]),c(0.025,0.5,0.975)))

smoothTraj<-function(line){
  name<-names(line)
  line<-as.numeric(line)
  first<-(line[1]+line[2])/2;last<-(line[length(line)-1]+line[length(line)])/2
  ret<-c(first,unlist(lapply(2:(length(line)-1),function(step){return((line[step-1]+line[step]+line[step+1])/3)})),
         last)
  names(ret)<-name;return(ret)}

ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")
Toplot<-Turn
windows()
plot(colnames(Toplot),Toplot[1,,"0.5"],type="n",ylab="",xlab="",ylim=c(min(Toplot),max(Toplot)))
invisible(lapply(1:length(treatments),function(tr){
    toplot<-Toplot[which(rownames(Toplot)%in%treatments[[tr]]),,]
    lapply(1:nrow(toplot),function(li){
      lines(colnames(toplot),smoothTraj(toplot[li,,"0.5"]),col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(toplot),rev(colnames(toplot))),
              c(smoothTraj(toplot[li,,"0.025"]),rev(smoothTraj(toplot[li,,"0.975"]))),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
}))
  
save(Nest,
     file="P:/Private/Taff/These/Redaction/DemographicDecomposition/RecruitmentTrajectories/DB/Nestedness_toInit")
save(Turn,
     file="P:/Private/Taff/These/Redaction/DemographicDecomposition/RecruitmentTrajectories/DB/Turnover_toInit")
save(Nest,file="Recruitment/DB/Nestedness_toInit")
save(Turn,file="Recruitment/DB/Turnover_toInit")


  
