setwd("P:/Private/Taff/These/Codes/DemographicDynamics")
load("DB/Paracou_R_Subdivided_ok")
source("Vernacular_handle.R")
source("TraitsMiceFilling.R")
library(cluster)
load("DB/Alpha_Plots")

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("Hmax","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")

#InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
#  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]return(ret[which(!duplicated(ret)),])}))
#InventorySp<-InventorySp[which(!duplicated(InventorySp)),]
#save(InventorySp,file="DB/InventorySp")
load(InventorySp,file="DB/InventorySp")

traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
traits_filled<-aggregate(traits_filled[,TraitsName],list(traits_filled$name),median)
rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]

#Plots of the different treatments
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

##Separate by year
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

Nrep<-50
RecPun_Fun<-lapply(1:Nrep,function(rep){
  Ret<-do.call(cbind,lapply(Recruits3,function(yr){
  ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
  ret<-yr[which(yr[,"n_parcelle"]==plot),]
  ret<-Replacement(ret,Alpha=alphas_plot[[which(names(alphas_plot)==plot)]])
  
  #traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
  #traits_filled<-aggregate(traits_filled[,TraitsName],list(traits_filled$name),median)
  #rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]
  tra<-traits_filled[which(rownames(traits_filled)%in%ret),];tra<-tra[order(rownames(tra)),]
  ret<-ret[which(ret%in%rownames(tra))]
  
  dissim<-as.matrix(daisy(tra,metric="gower"))
  dissim <- 1 - dissim/max(dissim)
  
  return(expq(Hqz(as.AbdVector(tapply(ret,ret,length)), q=2, dissim,Correction="None"),q=2))
  }))
  ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
  ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
  ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
  return(ret["ret"])
  }))
colnames(Ret)<-names(Recruits3)
return(Ret)})
RecPun_Fun<-array(unlist(RecPun_Fun),dim=c(nrow(RecPun_Fun[[1]]),ncol(RecPun_Fun[[1]]),Nrep),
           dimnames=list(rownames(RecPun_Fun[[1]]),colnames(RecPun_Fun[[1]]),1:Nrep))
RecPun_Fun<-lapply(c(0.025,0.5,0.975),function(quant){
  return(apply(RecPun_Fun,c(1,2),function(rep){return(quantile(rep,probs = quant,na.rm=T))}))
  })
RecPun_Fun<-array(unlist(RecPun_Fun),dim=c(nrow(RecPun_Fun[[1]]),ncol(RecPun_Fun[[1]]),3),
           dimnames=list(rownames(RecPun_Fun[[1]]),as.numeric(colnames(RecPun_Fun[[1]]))-1984,c(0.025,0.5,0.975)))

