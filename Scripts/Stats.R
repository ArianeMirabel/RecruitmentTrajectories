load("DB/Turnover_toInit")
treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

Maxs<-cbind(apply(Turn,1,max),treats[,"treat"])
colnames(Maxs)<-c("Max","treat")

cor.test(Maxs[,"Max"],Maxs[,"treat"],method="spearman")

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

c<-c("L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD","moisture","Bark_thick")

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

## Traits percent uncertainty
Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits1<-Traits1[which(apply(Traits1,1,function(li){return(any(!is.na(li)))})),]
Traits1["name"]<-as.factor(paste(Traits1$Genus,"_",Traits1$species,sep=""))
traits<-Traits1[,c("Family","Genus","species","name","bar_code","thickness","SPAD","toughness","dry_mass","traits_surf_area","ind_surf_area",
                    "sapwood_dens","moisture","bark_thick")]
colnames(traits)<-c("Family","Genus","species","name","bar_code","L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD",
                     "moisture","Bark_thick")
traits[which(!is.na(traits[,"LA"])),"SLA"]<-
  traits[which(!is.na(traits[,"LA"])),"LA"]/traits[which(!is.na(traits[,"LA"])),"L_DryMass"]
Seltraits<-c("L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")#"S_mass","moisture")
traits<-traits[,c("Family","Genus","name","bar_code",Seltraits)]

library(doBy)
traitsSp <- summaryBy(L_thickness + L_chloro + L_toughness + L_DryMass + SLA + WD + Bark_thick ~ name, 
                      data=traits, FUN = function(x){return(mean(x,na.rm=T))},keep.names = T)
length(which(apply(traitsSp[,Seltraits],1,anyNA)))/nrow(traitsSp)

length(which(apply(traits[,Seltraits],1,anyNA)))/nrow(traits)

load("DB/Paracou_R_Subdivided_ok")
diff <- setdiff(unique(Recruitment[,"name"]),unique(traits[,"name"]))
diff <- diff[which(!grepl("Indet.",diff))]
length(diff)/length(unique(Recruitment[,"name"]))
length(which(Recruitment[,"name"] %in% diff))/nrow(Recruitment)

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
load("DB/RecruitmentPunctual_Nullmodel_Diff");load("DB/RecruitmentPunctual_Functional_Nullmodel_Diff")

treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

Tax <- RecPun
Fun <- RecPun_Fun_DiffNull
  
for(q in 1:3){
  Ret<-Tax[[q]][,,"0.5"]   
  Ret <- apply(Ret,2,function(col){return(abs(col-Ret[,1]))})
  Ret<-apply(Ret,1,max)
  Ret<-apply(cbind(Ret,names(Ret)),2,as.numeric)
  colnames(Ret)<-c("Max","plot")
  Ret<-merge(Ret,AGBloss_cor,by="plot")
  assign(c("Richness","Shannon","Simpson")[q],Ret)
}

cor.test(Richness[,"Max"],Richness[,"AGB"],method="spearman")
cor(Shannon[,"Max"],Shannon[,"AGB"],method="spearman")
cor.test(Simpson[,"Max"],Simpson[,"AGB"],method="spearman")

Rao <- apply(Fun[,,"0.5"],2,function(col){return(abs(col- Fun[,1,"0.5"]))})
Rao <- apply(Rao,1,max)
Rao <- apply(cbind(Rao,names(Rao)),2,as.numeric)
colnames(Rao)<-c("Max","plot")
Rao<-merge(Rao,AGBloss_cor,by="plot")
cor.test(Rao[,"Max"],Rao[,"AGB"],method="spearman")


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

## Dominant species
Rec_stat <- Recruitment[,c("n_parcelle","name","campagne")]

Rec_statD <- subset(Rec_stat, n_parcelle %in% c(2,7,9,3,5,10,4,8,12))
Rec_statD1 <- Rec_stat[which(Rec_stat[,"campagne"]-1984 <=15),]
Dom_D1 <- sort(tapply(Rec_statD1$name,Rec_statD1$name,length),decreasing=T)[1:10]
Rec_statD2 <- Rec_stat[which(Rec_stat[,"campagne"]-1984 > 15),]
Dom_D2 <- sort(tapply(Rec_statD2$name,Rec_statD2$name,length),decreasing=T)[1:10]

length(grep("sp.",Recruitment$Espece))/nrow(Recruitment)*100

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
ColorsTr<-c("darkolivegreen2","gold","orangered","darkred")


Ntrait<-names(Rec_CWM)

windows()
par(mfrow=c(4,2))

PlotCWM<-function(TrajTraits){
  for (Ntrait in 1:length(TrajTraits)){
    toplot<-Rec_CWM[[Ntrait]][,,"0.5"]
    nm<-colnames(toplot)
    toplot<-as.data.frame(cbind(rownames(toplot),toplot))
    
    toplot<-reshape(toplot, direction="long", varying = 2:14, idvar="V1", v.names="val", timevar="time", times=nm)
    toplot<-as.data.frame(apply(toplot,2,as.numeric))
    plot(x=toplot[,"time"],y=toplot[,"val"],type="n",xlab="time",ylab=names(Rec_CWM)[Ntrait])
    
    lapply(1:4,function(treat){
      colTr<-col2rgb(ColorsTr[treat])
      colTr<-rgb(red=colTr[1], green=colTr[2], blue=colTr[3], maxColorValue = 255,alpha=75)
      dat<-toplot[which(toplot[,"V1"]%in%treatments[[treat]]),]
      points(x=dat[,"time"],y=dat[,"val"],col=colTr,pch=20,cex=0.5)
      gam_mod <- gam(val ~ s(time,k=5), data = dat, sp=0.1)
      xseq<-seq(from=5,to=30,length.out = 100)
      pred<-predict.gam(gam_mod,data.frame(time=xseq),type="response",se.fit=T)
      lines(y=pred$fit, x=xseq,col=ColorsTr[treat],lwd=1.5)
      polygon(c(xseq,rev(xseq)),c(pred$se.fit,rev(pred$se.fit)),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    })
  }
}


########## Breakpoint analysis
load("DB/RecruitmentPunctual");load("DB/RecruitmentPunctual_Nullmodel_Diff")
load("DB/RecruitmentPunctual_Functional");load("DB/RecruitmentPunctual_Functional_Nullmodel_Diff")
load("DB/Turnover_toInit")

load("DB/LostAGB")
norm<-(AGBloss_cor[,"AGB"]-min(AGBloss_cor[,"AGB"]))/(max(AGBloss_cor[,"AGB"])-min(AGBloss_cor[,"AGB"]))
ColorsDist <- colorRampPalette(c("darkolivegreen2","gold","orangered","darkred"))(12)[as.numeric(cut(norm, breaks = 12))]

RecPun_BP_3 <- do.call(rbind,lapply(1:length(RecPun), function(ind){
  
  Totest <- RecPun[[ind]][AGBloss_cor[,"plot"],,"0.5"]
  
  breakPoints <- do.call(cbind,lapply(1:nrow(Totest),function(li){
    totest <- Totest[li,]
    x <- as.numeric(names(totest))
    pairs <- combn(x,2)
    
    mse <- numeric(ncol(pairs))
    
    for(i in 1:ncol(pairs)){
      piecewise1 <- lm(totest ~ x*(x < pairs[1,i]) + x*(x < pairs[2,i] & x > pairs[1,i]) + x*( x >= pairs[2,i]))
      mse[i] <- summary(piecewise1)$sigma
    }
    mse <- as.numeric(mse)
    BPs <- pairs[,which(mse==min(mse))]
    
    piecewise2 <- lm(totest ~ x*(x <= BPs[1]) + x*(x < BPs[2] & x > BPs[1]) + x*( x >= BPs[2]))
    return(!is.na(piecewise2$coefficients["x:x >= pairs[2, i]TRUE"]))
  }))
  
  return(breakPoints)
}))
which(RecPun_BP_3)

RecPun_Fun_BP_3 <-unlist(lapply(1:nrow(RecPun_Fun),function(li){
    totest <- RecPun_Fun[AGBloss_cor[,"plot"],,"0.5"][li,]
    x <- as.numeric(names(totest))
    pairs <- combn(x,2)
    
    mse <- numeric(ncol(pairs))
    
    for(i in 1:ncol(pairs)){
      piecewise1 <- lm(totest ~ x*(x < pairs[1,i]) + x*(x < pairs[2,i] & x > pairs[1,i]) + x*( x >= pairs[2,i]))
      mse[i] <- summary(piecewise1)$sigma
    }
    mse <- as.numeric(mse)
    BPs <- pairs[,which(mse==min(mse))]
    
    piecewise2 <- lm(totest ~ x*(x <= BPs[1]) + x*(x < BPs[2] & x > BPs[1]) + x*( x >= BPs[2]))
    return(!is.na(piecewise2$coefficients["x:x >= pairs[2, i]TRUE"]))
  }))
which(RecPun_Fun_BP_3)

RecPun_BP_1 <- do.call(rbind,lapply(1:length(RecPun), function(ind){
  
  Totest <- RecPun[[ind]][AGBloss_cor[,"plot"],,"0.5"]
  
  breakPoints <- do.call(cbind,lapply(1:nrow(Totest),function(li){
    totest <- Totest[li,]
    x <- as.numeric(names(totest))
   
    mse <- numeric(length(x))
    
    for(i in 1:length(totest)){
      piecewise1 <- lm(totest ~ x*(x < x[i]) + x*(x>=x[i]))
      mse[i] <- summary(piecewise1)$sigma 
    }
    mse <- as.numeric(mse)
    mse <- which(mse==min(mse))
    
    wBP <- AIC(lm(totest ~ x*(x < x[mse]) + x*(x>=x[mse])))
    noBP <- AIC(lm(totest ~ x))

    return(wBP<noBP)
  }))
  
  return(breakPoints)
}))

windows()
layout(rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),c(0,0,4,4,0,0),c(0,0,4,4,0,0)))

RecPun_BP <- lapply(c(1,3), function(ind){
  
  Totest <- RecPun[[ind]][AGBloss_cor[,"plot"],,"0.5"]
  
  plot(colnames(Totest),Totest[1,], ylim=c(min(Totest),max(Totest)),type="n",ylab="diversity",xlab="year", main =names(RecPun)[ind])
  
  breakPoints <- do.call(cbind,lapply(1:nrow(Totest),function(li){
    totest <- Totest[li,]
    x <- as.numeric(names(totest))
    
    mse <- numeric(length(x))
    
    for(i in 1:length(x)){
      piecewise1 <- lm(totest ~ x*(x < x[i]) + x*(x>=x[i]))
      mse[i] <- summary(piecewise1)$sigma
    }
    mse <- as.numeric(mse)
    
    BP <- x[which(mse==min(mse))]
    
    piecewise2 <- lm(totest ~ x*(x < BP) + x*(x > BP))
    lmBP <- summary(piecewise2)$coefficients[,"Estimate"]
    Pval <- summary(piecewise2)$fstatistic
    Pval <- round(pf(Pval[1],Pval[2],Pval[3],lower.tail=F),3)
    
    lines(x,totest,col=ColorsDist[li])
    curve((lmBP["(Intercept)"] + lmBP["x < BPTRUE"]) + (lmBP["x"]+lmBP["x:x < BPTRUE"])*x, add=T, from= x[1], to=BP,col=ColorsDist[li])
    curve((lmBP["(Intercept)"] + lmBP["x > BPTRUE"]) + lmBP["x"]*x, add=T, from=BP, to=max(x),col=ColorsDist[li])
    abline(v=BP, lty=3,col=ColorsDist[li])
    return(c(BP,min(mse),Pval))
  }))
  
  rownames(breakPoints) <- c("BP","MSEmodel","Pval")  
  colnames(breakPoints) <- rownames(Totest)
  
  return(breakPoints)
})

plot(colnames(RecPun_Fun),RecPun_Fun[1,,"0.5"],
       ylim=c(min(RecPun_Fun[,,"0.5"]),max(RecPun_Fun[,,"0.5"])),type="n",ylab="diversity",xlab="year", main = "Rao")
RecPun_Fun_BP <- unlist(lapply(1:nrow(RecPun_Fun),function(li){
  
    totest <- RecPun_Fun[AGBloss_cor[,"plot"],,"0.5"][li,]
    x <- as.numeric(names(totest))
    
    mse <- numeric(length(x))
    
    for(i in 1:length(x)){
      piecewise1 <- lm(totest ~ x*(x < x[i]) + x*(x>=x[i]))
      mse[i] <- summary(piecewise1)$sigma
    }
    mse <- as.numeric(mse)
    
    BP <- x[which(mse==min(mse))]
    
    piecewise2 <- lm(totest ~ x*(x < BP) + x*(x > BP))
    lmBP <- summary(piecewise2)$coefficients[,"Estimate"]
    
    lines(x,totest,col=ColorsDist[li])
    curve((lmBP["(Intercept)"] + lmBP["x < BPTRUE"]) + (lmBP["x"]+lmBP["x:x < BPTRUE"])*x, add=T, from= x[1], to=BP,col=ColorsDist[li])
    curve((lmBP["(Intercept)"] + lmBP["x > BPTRUE"]) + lmBP["x"]*x, add=T, from=BP, to=max(x),col=ColorsDist[li])
    abline(v=BP, lty=3,col=ColorsDist[li])
    return(c(BP,min(mse)))
  }))


plot(colnames(Turn),Turn[1,,"0.5"],
     ylim=c(min(Turn[,,"0.5"]),max(Turn[,,"0.5"])),type="n",ylab="diversity",xlab="year", main = "Turnover")
RecPun_Turn_BP <- unlist(lapply(1:nrow(Turn),function(li){
  
  totest <- Turn[AGBloss_cor[,"plot"],,"0.5"][li,]
  x <- as.numeric(names(totest))
  
  mse <- numeric(length(x))
  
  for(i in 1:length(x)){
    piecewise1 <- lm(totest ~ x*(x < x[i]) + x*(x>=x[i]))
    mse[i] <- summary(piecewise1)$sigma
  }
  mse <- as.numeric(mse)
  
  BP <- x[which(mse==min(mse))]
  
  piecewise2 <- lm(totest ~ x*(x < BP) + x*(x > BP))
  lmBP <- summary(piecewise2)$coefficients[,"Estimate"]
  
  lines(x,totest,col=ColorsDist[li])
  curve((lmBP["(Intercept)"] + lmBP["x < BPTRUE"]) + (lmBP["x"]+lmBP["x:x < BPTRUE"])*x, add=T, from= x[1], to=BP,col=ColorsDist[li])
  curve((lmBP["(Intercept)"] + lmBP["x > BPTRUE"]) + lmBP["x"]*x, add=T, from=BP, to=max(x),col=ColorsDist[li])
  abline(v=BP, lty=3,col=ColorsDist[li])
  return(c(BP,min(mse)))
}))


### Uncertainty deviation
load("DB/RecruitmentPunctual");load("DB/RecruitmentPunctual_Nullmodel_Diff")

DevUncertainty_Taxo <-lapply(c(1,3),function(ind){
  recind<-RecPun[[ind]]
  ret<-mean(unlist(lapply(1:nrow(recind),function(Li){
    return(max(abs((recind[Li,,"0.975"]-recind[Li,,"0.025"])/recind[Li,,"0.5"]*100)))
  })))
  
})

DevUncertainty_Fun <-mean(unlist(lapply(1:nrow(RecPun_Fun),function(Li){
    return(max(abs((RecPun_Fun[Li,,"0.975"]-RecPun_Fun[Li,,"0.025"])/RecPun_Fun[Li,,"0.5"]*100)))
  })))

load("DB/Turnover_toInit")
DevUncertainty_Turn <-mean(unlist(lapply(1:nrow(Turn),function(Li){
  return(max(abs((Turn[Li,,"0.975"]-Turn[Li,,"0.025"])/Turn[Li,,"0.5"]*100)))
})))


