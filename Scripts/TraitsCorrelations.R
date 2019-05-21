library(lme4)
library(GeneNet)
library(metafor)

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("SLA","WD","Hmax","Bark_thick","L_toughness","L_chloro","L_thickness")#,"L_DryMass"

load("DB/InventorySp")

source("Scripts/TraitsMiceFilling.R")
  
Nrep<-100

Tcorr<-lapply(1:Nrep,function(rep){
  Traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
#Traits_filled<-aggregate(Traits_filled[,TraitsName],list(Traits_filled$name),median)
#rownames(Traits_filled)<-Traits_filled[,1]
traits_filled<-scale(Traits_filled[,TraitsName])
traits_filled<-as.data.frame(cbind(traits_filled,Traits_filled[,"name"]))
colnames(traits_filled)<-c(TraitsName,"name")

correl<-matrix(nrow=length(TraitsName),ncol=length(TraitsName),dimnames=list(TraitsName,TraitsName))
Pairs<-combn(TraitsName,2)
for (i in 1:ncol(Pairs)){
pair<-Pairs[,i]
fml<-as.formula(paste(pair[1], "~", pair[2], "+", "(1|name)"))
lmm<-lmer(fml, data=traits_filled)
correl[pair[1],pair[2]]<-fixef(lmm)[pair[2]]
}
return(z.transform(correl))
})

TraitCorr<-array(unlist(Tcorr),dim=c(nrow(Tcorr[[1]]),ncol(Tcorr[[1]]),Nrep),
      dimnames=list(TraitsName,TraitsName,1:Nrep))
TraitCorr<-apply(TraitCorr,c(1,2),function(c){return(round(transf.stor(mean(c,na.rm=T),2)))})
for (i in 1:length(TraitsName)){
  TraitCorr[TraitsName[i],TraitsName[i]]<-1
}
TraitCorr[which(is.na(TraitCorr))]<-""
as.data.frame(TraitCorr)

save(TraitCorr,file="DB/TraitsCorrelations")

Tcorr2<-lapply(1:Nrep,function(rep){
  Traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
  #Traits_filled<-aggregate(Traits_filled[,TraitsName],list(Traits_filled$name),median)
  #rownames(Traits_filled)<-Traits_filled[,1]
  traits_filled<-scale(Traits_filled[,TraitsName])
  traits_filled<-as.data.frame(cbind(traits_filled,Traits_filled[,"name"]))
  colnames(traits_filled)<-c(TraitsName,"name")
  
  ret<-corr.test(traits_filled[,TraitsName],use="pairwise.complete.obs",method="pearson")
  return(ret$r)
})

TraitCorr2<-array(unlist(Tcorr2),dim=c(nrow(Tcorr2[[1]]),ncol(Tcorr2[[1]]),Nrep),
                 dimnames=list(TraitsName,TraitsName,1:Nrep))
TraitCorr2<-apply(TraitCorr2,c(1,2),function(c){return(round(mean(c,na.rm=T),2))})


















