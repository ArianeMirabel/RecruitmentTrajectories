load("DB/Paracou_R_Subdivided_ok")
source("Scripts/Vernacular_handle.R")
source("Scripts/TraitsMiceFilling.R")
library(cluster)
load("DB/Alpha_Plots")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("SLA","WD","Hmax","Bark_thick","L_chloro","L_thickness")#,"L_DryMass","L_toughness"

#InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
#  ret<-ret[which(!grepl("Indet.",ret[,"name"])),];return(ret[which(!duplicated(ret)),])}))
#InventorySp<-InventorySp[which(!duplicated(InventorySp)),]
#save(InventorySp,file="DB/InventorySp")
load("DB/InventorySp")

#Plots of the different treatments
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

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

Nrep<-10
RecPun_Fun<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
  ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
  ret<-yr[which(yr[,"n_parcelle"]==plot),]
  ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
  #ret<-merge(as.data.frame(ret),InventorySp,by.x="ret",by.y="name",all.x=T)[,"Genre"]
  #ret<-as.character(ret[which(!is.na(ret))])
  
  #Traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
  #Traits_filled<-aggregate(Traits_filled[,TraitsName],list(Traits_filled$name),median)
  #rownames(Traits_filled)<-Traits_filled[,1];Traits_filled<-Traits_filled[,TraitsName]
  tra<-Traits_filled[which(rownames(Traits_filled)%in%ret),];tra<-tra[order(rownames(tra)),]
  ret<-ret[which(ret%in%rownames(tra))]
  
  dissim<-as.matrix(daisy(tra,metric="gower"))
  dissim <- 1 - dissim/max(dissim)
   
  return(expq(Hqz(as.AbdVector(tapply(ret,ret,length)), q=2, dissim,Correction="Best"),q=2))
  }))
  ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
  ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
  ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
  return(ret["ret"])
  }))
return(Ret)})
RecPun_Fun<-lapply(RecPun_Fun,function(rep){return(smooth(rep,2))}) # moving average, path=2
RecPun_Fun<-array(unlist(RecPun_Fun),dim=c(nrow(RecPun_Fun[[1]]),ncol(RecPun_Fun[[1]]),Nrep),
           dimnames=list(rownames(RecPun_Fun[[1]]),names(Recruits3),1:Nrep))
RecPun_Fun<-lapply(c(0.025,0.5,0.975),function(quant){
  return(apply(RecPun_Fun,c(1,2),function(rep){return(quantile(rep,probs = quant,na.rm=T))}))
  })
RecPun_Fun<-array(unlist(RecPun_Fun),dim=c(nrow(RecPun_Fun[[1]]),ncol(RecPun_Fun[[1]]),3),
           dimnames=list(rownames(RecPun_Fun[[1]]),as.numeric(colnames(RecPun_Fun[[1]]))-1984,c(0.025,0.5,0.975)))

save(RecPun_Fun,file="DB/RecruitmentPunctual_Functional")

#Missings<-Recruitment[which(Recruitment[,"name"]%in%missings),]
#Missings<-Missings[which(Missings$n_parcelle%in%1:12),]
#Missings<-Missings[which(!Missings$campagne%in%c(1998,2000,2002,2004,2006)),]

#lapply(unique(Missings[,"campagne"]),function(yr){Yr<-Missings[which(Missings$campagne==yr),];return(tapply(Yr$n_parcelle,Yr$n_parcelle,length))})

#plot(1:12,tapply(Missings$n_parcelle,Missings$n_parcelle,length),type="n")
#lapply(1:4,function(tr){
#  toplot<-Missings[which(Missings$n_parcelle%in%treatments[[tr]]),]
#  points(treatments[[tr]],tapply(toplot$n_parcelle,toplot$n_parcelle,length),
#         col=c("darkolivegreen2","deepskyblue2","darkorange1","red2")[tr],pch=19)
#})

#toplot<-tapply(Missings$campagne,Missings$campagne,length)
#plot(names(toplot),toplot)


###########################################
RecPun_Fun_Null<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
    ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
      ret<-yr[which(yr[,"n_parcelle"]==plot),]
      ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
      #ret<-merge(as.data.frame(ret),InventorySp,by.x="ret",by.y="name",all.x=T)[,"Genre"]
      #ret<-as.character(ret[which(!is.na(ret))])
      
      Traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
      Traits_filled<-aggregate(Traits_filled[,TraitsName],list(Traits_filled$name),median)
      rownames(Traits_filled)<-Traits_filled[,1];Traits_filled<-Traits_filled[,TraitsName]
      tra<-Traits_filled;rownames(tra)<-rownames(tra)[sample(nrow(tra))]
      tra<-Traits_filled[which(rownames(Traits_filled)%in%ret),];tra<-tra[order(rownames(tra)),]
      ret<-ret[which(ret%in%rownames(tra))]
      
      dissim<-as.matrix(daisy(tra,metric="gower"))
      dissim <- 1 - dissim/max(dissim)
      
      return(expq(Hqz(as.AbdVector(tapply(ret,ret,length)), q=2, dissim,Correction="Best"),q=2))
    }))
    ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
    ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
    ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
    return(ret["ret"])
  }))
  colnames(Ret)<-names(Recruits3)
  return(Ret)})
RecPun_Fun_Null<-lapply(RecPun_Fun_Null,function(rep){return(smooth(rep,2))}) # moving average, path=2
RecPun_Fun_Null<-array(unlist(RecPun_Fun_Null),dim=c(nrow(RecPun_Fun_Null[[1]]),ncol(RecPun_Fun_Null[[1]]),Nrep),
                  dimnames=list(rownames(RecPun_Fun_Null[[1]]),colnames(RecPun_Fun_Null[[1]]),1:Nrep))
RecPun_Fun_Null<-lapply(c(0.025,0.5,0.975),function(quant){
  return(apply(RecPun_Fun_Null,c(1,2),function(rep){return(quantile(rep,probs = quant,na.rm=T))}))
})
RecPun_Fun_Null<-array(unlist(RecPun_Fun_Null),dim=c(nrow(RecPun_Fun_Null[[1]]),ncol(RecPun_Fun_Null[[1]]),3),
                  dimnames=list(rownames(RecPun_Fun_Null[[1]]),as.numeric(colnames(RecPun_Fun_Null[[1]]))-1984,c(0.025,0.5,0.975)))

save(RecPun_Fun_Null,file="DB/RecruitmentPunctual_Functional_Nullmodel")
