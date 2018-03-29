setwd("P:/Private/Taff/These/Codes/DemographicDynamics")
load("DB/Paracou_R_Subdivided_ok")
library(RColorBrewer)
library(reshape2)
library(cluster)

source("Vernacular_handle.R")
source("Overview_functions.R")
source("TraitsMiceFilling.R")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("Hmax","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")

load("DB/Alpha_Plots")
InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

#Plots of the different treatments
T0<-c(1,6,11); T1<-c(2,7,9); T2<-c(3,5,10); T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

##Separate by year
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
Traits_traj<-lapply(TraitsName,function(traitName){
  Traj<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
    ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
      ret<-yr[which(yr[,"n_parcelle"]==plot),]
      ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
      
      #Traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
      traits_filled<-aggregate(Traits_filled[,traitName],list(Traits_filled$name),median)
      rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled["x"]
      
      tra<-traits_filled[which(rownames(traits_filled)%in%ret),,drop=F];tra<-tra[order(rownames(tra)),,drop=F]
      ret<-ret[which(ret%in%rownames(tra))]
      ret<-as.ProbaVector(tapply(ret,ret,length))
      
      return(sum(tra*ret))
    }))
    ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
    ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
    ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
    return(ret["ret"])
  }))
  return(Ret)})
Traj<-lapply(Traj,function(rep){return(smooth(rep,2))}) # moving average, path=2
Traj<-array(unlist(Traj),dim=c(nrow(Traj[[1]]),ncol(Traj[[1]]),Nrep),
                  dimnames=list(rownames(Traj[[1]]),names(Recruits3),1:Nrep))
Traj<-lapply(c(0.025,0.5,0.975),function(quant){
  return(apply(Traj,c(1,2),function(rep){return(quantile(rep,probs = quant,na.rm=T))}))
})
Traj<-array(unlist(Traj),dim=c(nrow(Traj[[1]]),ncol(Traj[[1]]),3),
                  dimnames=list(rownames(Traj[[1]]),as.numeric(colnames(Traj[[1]]))-1984,c(0.025,0.5,0.975)))
})
names(Traits_traj)<-TraitsName

save(Traits_traj,file="DB/Recruits_CWM_sansHmax")
save(Rec_CWM,file="P:/Private/Taff/These/Redaction/3_RecruitmentTrajectories/DB/Recruits_CWM")
load("DB/Recruits_CWM_Hmax")

ColorsTr<-c("chartreuse3","deepskyblue2","darkorange1","red2")

windows()
par(mfrow=c(2,4),mar=c(2,2,2,1))
for (Ntrait in names(Traits_traj)){
  Toplot<-Traits_traj[[Ntrait]]
  plot(colnames(Toplot),Toplot[1,,"0.5"],type="n",xaxt="n",xlab="",
       ylab="",cex.lab=1.5,ylim=c(min(Toplot),max(Toplot)))
  axis(1,at=seq(5,30,5),labels=T) 
  mtext(Ntrait,3,adj=0,line = 0.5,cex=0.8)

invisible(lapply(1:length(treatments),function(tr){
  toplot<-Toplot[which(rownames(Toplot)%in%treatments[[tr]]),,]  
  lapply(1:nrow(toplot),function(Li){
   lines(colnames(toplot),toplot[Li,,"0.5"],col=ColorsTr[[tr]],lwd=2)
   polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[Li,,"0.025"],rev(toplot[Li,,"0.975"])),
            col=rgb(0,0,0,alpha=0.05),border=NA)})
}))
}

windows()
par(mfrow=c(2,4),mar=c(2,2,2,1))
for (Ntrait in names(Traits_traj)){
  Toplot<-Traits_traj[[Ntrait]]
  plot(colnames(Toplot),Toplot[1,,"0.5"],type="n",xaxt="n",xlab="",
       ylab="",cex.lab=1.5,ylim=c(min(Toplot),max(Toplot)))
  axis(1,at=seq(5,30,5),labels=T) 
  mtext(Ntrait,3,adj=0,line = 0.5,cex=0.8)
  
  invisible(lapply(1:length(treatments),function(tr){
    toplot<-Toplot[which(rownames(Toplot)%in%treatments[[tr]]),,]  
      lines(colnames(toplot),apply(toplot[,,"0.5"],2,median),col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(toplot),rev(colnames(toplot))),
              c(apply(toplot[,,"0.025"],2,median),rev(apply(toplot[,,"0.975"],2,median))),
              col=rgb(0,0,0,alpha=0.05),border=NA)
  }))
}

  