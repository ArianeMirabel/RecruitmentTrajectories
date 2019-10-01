load("DB/Turnover_toInit")
treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

Maxs<-cbind(apply(Turn,1,max),treats[,"treat"])
colnames(Maxs)<-c("Max","treat")

cor(Maxs[,"Max"],Maxs[,"treat"],method="spearman")

##################################################"

Traits<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits["name"]<-as.factor(paste(Traits$Genus,"_",Traits$species,sep=""))
traits<-Traits1[,c("Family","Genus","species","name","bar_code","thickness","SPAD","toughness","dry_mass","traits_surf_area","ind_surf_area",
                   "sapwood_dens","moisture","bark_thick")]
colnames(traits)<-c("Family","Genus","species","name","bar_code","L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD",
                    "moisture","Bark_thick")
##SLA calc for empty values with LA available
traits[which(!is.na(traits[,"LA"])),"SLA"]<-
  traits[which(!is.na(traits[,"LA"])),"LA"]/traits[which(!is.na(traits[,"LA"])),"L_DryMass"]

traitsName<-c("L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD","moisture","Bark_thick")

allSd<-apply(traits[,traitsName],2,function(x){return(sd(x, na.rm=T))})

SpSd<-aggregate(traits[,traitsName],list(traits$name),function(x){return(sd(x, na.rm=T))})
SpSd<-apply(SpSd[,traitsName],2,function(x){return(sd(x, na.rm=T))})

GenSd<-aggregate(traits[,traitsName],list(traits$Genus),function(x){return(sd(x, na.rm=T))})
GenSd<-apply(GenSd[,traitsName],2,function(x){return(sd(x, na.rm=T))})

Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
Traits2<-Traits2[which(apply(Traits2,1,function(li){return(any(!is.na(li)))})),]
Traits2["name"]<-sub(" ","_",Traits2[,"Name"])
traits2<-Traits2[,c("name","Hauteur","Masse")]
colnames(traits2)<-c("name","Hmax","S_mass")
traits2[which(traits2[,"name"]=="Sterculia_speciosa"),"Hmax"]<-43

traits<-merge(traits2,traits,by="name")[,c("Hmax","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")]
### To be continued


#########################################"
# Test on traits database

load("DB/BridgeData")

Phylo<-"Family" #"Genus" 

# Par famille ou par genre, comparer la variance des communautés à une distribution normale
Sep_t<-lapply(c("L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Hmax"),
              function(trait){
  sep<-lapply(as.character(unique(traits[,Phylo])),function(Fam){
    ret<-traits[which(traits[,Phylo]==Fam),trait]
    if(length(na.omit(ret))>1){
    ret<-var.test(na.omit(ret),na.omit(traits[,trait]))
    return(c(ret$p.value,ret$estimate,Fam))}
    })
  sep<-do.call(rbind,sep)
  colnames(sep)<-c("Pval","ratio",'phylo')
  return(sep)
})
names(Sep_t)<-c("L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Hmax")
lapply(Sep_t,function(sep){return(sep[which(as.numeric(sep[,2])>1),])})
  
# Correlations entre traits
Sp_means<-lapply(c("L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Hmax"), 
                 function (Tr){
Sp_mean<-by(data=traits[,Tr],traits[,"name"],function(x){return(mean(x, na.rm=T))},simplify=F)
return(unlist(Sp_mean))
})
Sp_means<-do.call(cbind,Sp_means)
colnames(Sp_means)<-c("L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Hmax")
Sp_means<-Sp_means[which(!apply(Sp_means,1,anyNA)),]

Corr<-cor(Sp_means)
Corr<-round(Corr,2)
save(Corr,file="DB/Traitscorrelations")
Corr<-cbind(rownames(Corr),Corr)
write.table(Corr,file="DB/Traitscorrelations.csv",sep=";",row.names=F)
  
############################################
## Spearman tests

load("DB/RecruitmentPunctual");load("DB/RecruitmentPunctual_Functional");load("DB/LostAGB")

treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

for(q in 1:3){
  Ret<-RecPun[[q]][,,"0.5"]     
  Ret<-apply(Ret,1,max)
  Ret<-apply(cbind(Ret,names(Ret)),2,as.numeric)
  colnames(Ret)<-c("Max","plot")
  Ret<-merge(Ret,AGBloss,by="plot")
  assign(c("Richness","Shannon","Simpson")[q],Ret)
}

cor(Richness[,"Max"],Richness[,"AGB"],method="spearman")
cor(Shannon[,"Max"],Shannon[,"AGB"],method="spearman")
cor(Simpson[,"Max"],Simpson[,"AGB"],method="spearman")

cor(Simpson[which(!Simpson$plot%in%c(8,12)),"Max"],Simpson[which(!Simpson$plot%in%c(8,12)),"treat"],method="spearman")

## Species repartition

load("DB/Paracou_R_Subdivided_ok")
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

Rec_stat<-Recruitment[,c("n_parcelle","name")]
Rec_stat<-Rec_stat[which(!duplicated(Rec_stat)),]
Rec_stat2<-unlist(lapply(unstack(Rec_stat),function(x){
  return(length(which(x%in%1:12)))}))

length(which(Rec_stat2==1))
       
Rec_stat3<-lapply(treatments,function(tr){
  ret<-Rec_stat[which(Rec_stat[,"n_parcelle"]%in%tr),]
  ret<-unstack(ret)[which(!unlist(lapply(unstack(ret),is.null)))]
  #return(length(ret))
  return(names(ret))
})

#Species only recruited in control plots
Only0<-setdiff(setdiff(setdiff(Rec_stat3[[1]],Rec_stat3[[2]]),Rec_stat3[[3]]),Rec_stat3[[4]])
#Species only recruited in T1
Only1<-setdiff(setdiff(setdiff(Rec_stat3[[2]],Rec_stat3[[1]]),Rec_stat3[[3]]),Rec_stat3[[4]])
#Species only recruited in T2
Only2<-setdiff(setdiff(setdiff(Rec_stat3[[3]],Rec_stat3[[1]]),Rec_stat3[[2]]),Rec_stat3[[4]])
#Species only recruited in T3
Only3<-setdiff(setdiff(setdiff(Rec_stat3[[4]],Rec_stat3[[1]]),Rec_stat3[[3]]),Rec_stat3[[3]])

# NMDS of recruitment
library(ade4)
library(factoextra)

ColorsTr<-c("darkolivegreen2","gold","orangered","darkred")

dates<-sort(names(LivingStand_all))[-1]
dates<-dates[which(!dates%in%c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016"))]

all<-as.character(Recruitment[,"name"])
all<-as.data.frame(as.ProbaVector(tapply(all,all,length)))

Mat<-lapply(1:12,function(p){
    temp<-Recruitment[which(Recruitment[,"n_parcelle"]==p),]
      if(nrow(temp)==0){return(NA)}
      if(nrow(temp)!=0){temp<-as.character(temp[!duplicated(temp),"name"])
      return(as.ProbaVector(tapply(temp,temp,length)))}
})
names(Mat)<-1:12

DistMat<-do.call(cbind,lapply(Mat,function(plo){
  plo<-as.data.frame(plo)
  ret<-merge(all,plo,by="row.names",all.x=T)
  ret[is.na(ret)]<-0
  rownames(ret)<-ret[,"Row.names"];ret<-ret["plo"]
  return(ret)#dist(t(ret[,c(2:3)]))
  }))
colnames(DistMat)<-1:12

AFC<-dudi.coa(DistMat, scannf=F,nf=2)
plot(AFC)

Pcoord<-AFC$co;rownames(Pcoord)<-1:12
Spcont<-get_ca_row(AFC)$contrib
Sp_1<-rownames(Spcont[order(Spcont[,"Dim.1"],decreasing=T)[1:10],])
Sp_2<-rownames(Spcont[order(Spcont[,"Dim.2"],decreasing=T)[1:10],])
SpCont<-AFC$li[unique(c(Sp_1,Sp_2)),]

windows()
plot(Pcoord,type="n")
lapply(1:4,function(tr){
 Tr<-treatments[[tr]]
 text(x=Pcoord[Tr,1],y=Pcoord[Tr,2],col=ColorsTr[tr],labels=rownames(Pcoord[Tr,]),cex=1.5)
})
text(x=SpCont[,"Axis1"],y=SpCont[,"Axis2"],labels=rownames(SpCont),cex=0.8)


############## Number of recruits
load("DB/Paracou_R_Subdivided_ok")

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

RecCount<-do.call(cbind,lapply(Recruits3,function(yr){
      ret<-unlist(lapply(sort(unique(yr[,"n_parcelle"])),function(plot){
        return(nrow(yr[which(yr[,"n_parcelle"]==plot),]))
      }))
      ret<-as.data.frame(ret,row.names=as.character(sort(unique(yr[,"n_parcelle"]))))
      ret<-merge(ret,as.data.frame(1:12,row.names=as.character(1:12)),by="row.names",all.y=TRUE)[,1:2]
      ret<-ret[order(as.numeric(ret[,"Row.names"])),];rownames(ret)<-ret[,"Row.names"]
      return(ret["ret"])
    }))
RecCount<-smooth(recind,2) # moving average, path=2
colnames(RecCount)<-as.numeric(names(Recruits3))-1984

windows()
ColorsTr<-c("darkolivegreen2","gold","orangered","darkred")
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12);treatments<-list(T0,T1,T2,T3)

Ylim=c(min(RecCount),max(RecCount))
plot(colnames(RecCount),RecCount[1,],type="n",ylim=Ylim,xlab="",ylab="",xaxt="n")
    invisible(lapply(1:length(treatments),function(tr){
      toplot<-RecCount[which(rownames(RecCount)%in%treatments[[tr]]),]
      lapply(1:nrow(toplot),function(Li){
        li<-toplot[Li,]
        lines(colnames(toplot),toplot[Li,],col=ColorsTr[[tr]],lwd=2)
        points(x=names(li[which(li==max(li))]),y=max(li),pch=16,col=ColorsTr[[tr]])
        #text(labels=max(li),x=5,y=max(li),col=ColorsTr[[tr]])
        lines(x=c(3,names(li[which(li==max(li))])),y=rep(max(li),2),col=ColorsTr[[tr]],lty=3)
})
        boxed.labels(x=5,y=max(toplot),labels=max(toplot),col=ColorsTr[[tr]],bg="white",border=F)
}))
axis(1,at=seq(5,30,5),labels=T)

MeanEff<-unlist(lapply(treatments,function(tr){
  mean(RecCount[which(rownames(RecCount)%in%tr),])}))
names(MeanEff)<-names(treatments)


#### GAM for CWM trajectories
library("mgcv")
load("DB/Recruits_CWM")

T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3);names(treatments)<-c("0","1","2","3")

Ntrait<-names(Rec_CWM)

windows()
par(mfrow=c(4,2))
GAM_T0<-lapply(1:length(Rec_CWM),function(tra){
  treat<-lapply(1:4,function(tre){
  toplot<-Rec_CWM[[tra]][treatments[[1]],,"0.5"]
  toplot<-as.data.frame(cbind(rownames(toplot),toplot))
  toplot<-reshape(toplot,direction="long",idvar="rownames(toplot)",varying=list(2:ncol(toplot)))[,2:3]
  toplot<-as.data.frame(apply(toplot,2,as.numeric))
  colnames(toplot)<-c("time","div")
  gam_mod <- gam(div ~ s(time,k=3), data = toplot, sp=0.1)
  return(gam_mod)    
  })
  lapply(1:4,function(tre){
    plot(treat[[tre]], col=ColorsTr[tre],residuals=F)
    pred<-predict(treat[[3]], type = "terms",se.fit=T)
    xseq<-seq(1,10,length.out = length(pred$fit))
    lines(y=pred$fit,x=xseq)
    polygon(c(xseq,rev(xseq)),c(pred$se.fit,rev(pred$se.fit)),
            col=rgb(0,0,0,alpha=0.05),border=NA)
  
    par(new=F)
    plot(treat[[2]], col=ColorsTr[tre],residuals=F)
    par(new=F)
    plot(treat[[3]], col=ColorsTr[tre],residuals=F)
    
    
  })
  par(new=T)
return(treat)
})

ColorsTr<-c("darkolivegreen2","gold","orangered","darkred")

lapply(1:4,function(tre){
  plot(treat[[tre]], col=ColorsTr[tre])
  })

plot(toplot)
gam_mod <- gam(div ~ s(time,k=3), data = toplot, sp=0.1)

# Plot the results
plot(gam_mod, residuals = TRUE, pch = 1)




