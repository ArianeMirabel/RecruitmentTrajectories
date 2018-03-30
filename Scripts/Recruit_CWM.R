load("DB/Paracou_R_Subdivided_ok")
library(reshape2)
library(cluster)

source("Scripts/Vernacular_handle.R")
source("Scripts/Overview_functions.R")
source("Scripts/TraitsMiceFilling.R")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("L_thickness","L_chloro","L_toughness","SLA","WD","Bark_thick","Hmax")#"L_DryMass",

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
Rec_CWM<-lapply(TraitsName,function(traitName){
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
names(Rec_CWM)<-TraitsName

save(Rec_CWM,file="P:/Private/Taff/These/Redaction/3_RecruitmentTrajectories/DB/Recruits_CWM")

