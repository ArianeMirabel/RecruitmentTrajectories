library(ade4)
library(entropart)
library(vegan)

load("DB/Paracou_R_Subdivided_ok")

load("DB/BotaVern")
load("DB/Alpha_Plots")
source("Scripts/Vernacular_handle.R")

Verlist<-c("muamba","win udu","aganiamai","buguni","patawa")

AlphaRec<-alpha_construct(Life)
Nrep<-50
trials<-unlist(lapply(Verlist,function(name){
  ret<-unlist(lapply(1:Nrep,function(rep){
   return(Dirichlet_draw(AlphaRec[,name]))}))
  ret<-tapply(ret,ret,length)
  return(names(ret[which(ret==max(ret))]))
  }))
  



T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

dates<-unique(Recruitment[,"campagne"])[seq(2,29,2)]

AllInv<-unique(Recruitment[,"name"])
AllInv<-as.data.frame(AllInv[which(AllInv!="Indet._Indet.")])
rownames(AllInv)<-AllInv[,1]

RecList<-lapply(1:12,function(p){
    Traj<-Recruitment[which(Recruitment[,"n_parcelle"]==p),]
    Traj<-Traj[which(Traj[,"name"]!="Indet._Indet."),"name"]
    Traj<-as.data.frame(tapply(Traj,Traj,length))
    Traj<-merge(AllInv,Traj,by="row.names",all.x=TRUE)
    Traj[which(is.na(Traj[,3])),3]<-0
    return(Traj[3])})
p<-1; Traj<-Recruitment[which(Recruitment[,"n_parcelle"]==p),]
Traj<-Traj[which(Traj[,"name"]!="Indet._Indet."),"name"]
Traj<-as.data.frame(tapply(Traj,Traj,length))
Traj<-merge(AllInv,Traj,by="row.names",all.x=TRUE)
Traj[which(is.na(Traj[,3])),3]<-0
Nm<-Traj[,1]
RecList<-do.call(cbind,RecList)
rownames(RecList)<-Nm
colnames(RecList)<-1:12
#RecList<-apply(RecList,2,as.ProbaVector)

RecList_last5<-lapply(1:12,function(p){
  Traj<-Recruitment[which(
    Recruitment[,"n_parcelle"]==p &
      Recruitment[,"campagne"]%in% c(2009,2011,2013,2015)),]
  Traj<-Traj[which(Traj[,"name"]!="Indet._Indet."),"name"]
  Traj<-as.data.frame(tapply(Traj,Traj,length))
  Traj<-merge(AllInv,Traj,by="row.names",all.x=TRUE)
  Traj[which(is.na(Traj[,3])),3]<-0
  return(Traj[3])})
RecList_last5<-do.call(cbind,RecList_last5)
rownames(RecList_last5)<-Nm
colnames(RecList_last5)<-1:12

Dom_last5<-lapply(treatments,function(tr){
  return(names(head(
    sort(rowSums(RecList_last5[,tr]),decreasing=T),n=5L)))
})
unique(unlist(Dom_last5[2:4]))
Dom_last5[[1]]

RecList_first5<-lapply(1:12,function(p){
  Traj<-Recruitment[which(
    Recruitment[,"n_parcelle"]==p &
      Recruitment[,"campagne"]%in% c(1994, 1997, 1999, 2001, 2003)),]
  Traj<-Traj[which(Traj[,"name"]!="Indet._Indet."),"name"]
  Traj<-as.data.frame(tapply(Traj,Traj,length))
  Traj<-merge(AllInv,Traj,by="row.names",all.x=TRUE)
  Traj[which(is.na(Traj[,3])),3]<-0
  return(Traj[3])})
RecList_first5<-do.call(cbind,RecList_first5)
rownames(RecList_first5)<-Nm
colnames(RecList_first5)<-1:12

Dom_first5<-lapply(treatments,function(tr){
  return(names(head(
    sort(rowSums(RecList_first5[,tr]),decreasing=T),n=5L)))
})
unique(unlist(Dom_first5[2:4]))
Dom_first5[[1]]

# Cluster de plot: clear separation only between control and disturbed
#Measure indicative value
Cont<-unlist(treatments[1])
Dist<-unlist(treatments[c(2,3,4)])

ControlA<-apply(RecList[,Cont],1,mean)/(apply(RecList[,Cont],1,mean)+apply(RecList[,Dist],1,mean))
ControlB<-apply(RecList[,Cont],1,function(row){
  length(which(row!=0))/length(row)})
IndVal_control<-round(ControlA*ControlB*100)
IndVal_control[is.na(IndVal_control)]<-0

DisturbA<-apply(RecList[,Dist],1,mean)/(apply(RecList[,Cont],1,mean)+apply(RecList[,Dist],1,mean))
DisturbB<-apply(RecList[,Dist],1,function(row){
  length(which(row!=0))/length(row)})
IndVal_disturb<-round(DisturbA*DisturbB*100)
IndVal_disturb[is.na(IndVal_disturb)]<-0

Index<-cbind(IndVal_control,IndVal_disturb);colnames(Index)<-c("Control",'Disturb')
Icontrol<-IndVal_control[which(IndVal_control>90)]
Idisturb<-IndVal_disturb[which(IndVal_disturb>90)]

# permutation of sites in site groups to test the significance
Nrep<-50
Permute<-lapply(1:Nrep,function(rep){
  Cont<-floor(runif(3,min=1,max=12))
  Dist<-setdiff(1:12,Cont)
  
  ControlA<-apply(RecList[,Cont],1,mean)/(apply(RecList[,Cont],1,mean)+apply(RecList[,Dist],1,mean))
  ControlB<-apply(RecList[,Cont],1,function(row){
    length(which(row!=0))/length(row)})
  IndVal_control<-round(ControlA*ControlB*100)
  IndVal_control[is.na(IndVal_control)]<-0
  
  DisturbA<-apply(RecList[,Dist],1,mean)/(apply(RecList[,Cont],1,mean)+apply(RecList[,Dist],1,mean))
  DisturbB<-apply(RecList[,Dist],1,function(row){
    length(which(row!=0))/length(row)})
  IndVal_disturb<-round(DisturbA*DisturbB*100)
  IndVal_disturb[is.na(IndVal_disturb)]<-0
  
  Index<-cbind(IndVal_control,IndVal_disturb);colnames(Index)<-c("Control",'Disturb')
  return(Index)
})
Permute<-array(unlist(Permute),dim=c(nrow(Permute[[1]]),ncol(Permute[[1]]),length(Permute)),
               dimnames=list(rownames(Permute[[1]]),colnames(Permute[[1]]),1:length(Permute)))
MuPerm<-apply(Permute,c(1,2),mean)
SdPerm<-apply(Permute,c(1,2),sd)

#Z-statistic
Zeta<- (Index - MuPerm) / SdPerm
PvalZ<-pnorm(-abs(Zeta))
# p-value, if between 
pnorm(-abs(Zeta[names(Icontrol),]))
pnorm(-abs(Zeta[names(Idisturb),]))

#T test
Ttest<-cbind(unlist(lapply(1:dim(Permute)[1],function(sp){
  return((t.test(Permute[sp,"Control",],mu=Index[sp,"Control"])$p.value))})),
  unlist(lapply(1:dim(Permute)[1],function(sp){
    return((t.test(Permute[sp,"Disturb",],mu=Index[sp,"Disturb"])$p.value))})))
rownames(Ttest)<-rownames(Permute)
colnames(Ttest)<-colnames(Permute)
PvalT_cont<-Ttest[names(Icontrol),]
PvalT_dist<-Ttest[names(Idisturb),]

indic<-Recruitment[which(Recruitment[,"name"]%in%names(Idisturb)),]

## NMDS of plots from indicative species inventory
Mat<-lapply(1:12,function(p){
  rec<-as.character(indic[which(indic[,"n_parcelle"]==p),"name"])
  rec<-as.data.frame(as.ProbaVector(tapply(rec,rec,length)))
  rec<-merge(data.frame(row.names = names(Idisturb)),rec,by="row.names",all=T)
  rec[is.na(rec)]<-0
  return(rec[,2])
  })
Mat<-dist(t(do.call(cbind,Mat)))
Mat<-metaMDS(Mat,distance="bray")

ColorsTr<-c("darkolivegreen2","gold","orangered","darkred")
plot(Mat,type="n")
lapply(1:length(treatments),function(tr){
  toplot<-Mat$points[treatments[[tr]],]
  text(toplot,labels=rownames(toplot),col=ColorsTr[tr])
})

## NMDS of years from indicativs species inventory
MatD<-lapply(sort(unique(indic[,"campagne"])),function(d){
  rec<-as.character(indic[which(indic[,"campagne"]==d),"name"])
  rec<-as.data.frame(tapply(rec,rec,length))
  #rec<-as.data.frame(as.ProbaVector(tapply(rec,rec,length)))
  rec<-merge(data.frame(row.names = names(Idisturb)),rec,by="row.names",all=T)
  rec[is.na(rec)]<-0
  return(rec[,2])
})
MatD<-do.call(cbind,MatD);colnames(MatD)<-sort(unique(indic[,"campagne"]));rownames(MatD)<-names(Idisturb)
#MatD<-dist(t(do.call(cbind,MatD)))
MatD_nmds<-metaMDS(t(MatD),distance="bray")

k<-0
while(!MatD_nmds$converged & k<=100){MatD_nmds<-invisible(metaMDS(MatD,distance="bray", previous.best = MatD_nmds));k<-k+1;print(k)}

ColorsD<-colorRampPalette(c("coral","cornflowerblue")) (length(unique(indic[,"campagne"])))
plot(MatD_nmds,type="n")
text(MatD_nmds$species,labels=rownames(MatD_nmds$species),col=ColorsD,cex=1.1)
text(MatD_nmds$points,labels=rownames(MatD_nmds$points),cex=0.9)
  
lapply(1:length(unique(indic[,"campagne"])),function(d){
  toplot<-Mat$points[treatments[[tr]],]
})

#Number of clusters: plot within group sum of square against number of clusters
maxClust<-8
wss <- (nrow(RecList)-1)*sum(apply(RecList,2,var))
for (i in 2:maxClust) wss[i] <- sum(kmeans(RecList,centers=i)$withinss)

plot(1:maxClust, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

# Statistical gap method, compare observed inv. and uniformely distributed frequencies
Nrep<-5
stat<-do.call(rbind,lapply(1:Nrep,function(rep){
  NullList<-do.call(cbind,lapply(1:12,function(p){
  return(runif(nrow(RecList),min=0,max=1))}))
  
  wssNull <- (nrow(NullList)-1)*sum(apply(NullList,2,var))
  for (i in 2:maxClust) wssNull[i] <- sum(kmeans(NullList,centers=i,iter.max=30)$withinss)
  
  return(log(wssNull)-log(wss))
}))
Gap<-colSums(stat)/Nrep
Std<-apply(stat,2,sd)

#The best number of gaps: gap(k) > gap(k+1)-std(k+1)
score <- 0
for (i in 2:length(Gap)-1) score[i] <-Gap[i+1]-Std[i+1]
score[length(Gap)]<-0
GapStat<-Gap-score

plot(1:14, GapStat[1:14], type="b", xlab="Nb of Clusters",
     ylab="Gap statistic")
abline(a=0,b=0,col='red')

Nclust<-which(GapStat>0)[1]

d<-dist(t(RecList)) # Euclidean distance between recruited communities based on inventory data
clusters<-hclust(d,"ward.D") 
kclust<-kmeans(d,centers=Nclust)
# Cluster de plot: clear separation only between control and disturbed
 
ACP<-dudi.pca(RecList,scannf = F,nf=3)

col<-c("darkolivegreen2","gold","orangered","darkred")

windows()
plot(ACP$co[,c("Comp1","Comp2")], type="n",xlim=c(-1,1),ylim=c(-1,1),xlab="",ylab="",
     pty="s",bty="n",axes=F,cex=1.3)
mtext("Axis  1",side=4,las=1,at=-0.05,line=-2);
mtext("Axis  2",side=3)
lapply(1:length(treatments),function(tr){
  Tr<-treatments[[tr]]
text(ACP$co[Tr,1], ACP$co[Tr,2], labels = rownames(ACP$co[Tr,]),col=col[tr],cex=1.6)
})
abline(h=0,v=0)

plot(ACP$co[,1], xlim=c(-1,-0.5),ylim=c(0,0), axes=FALSE, type = "n", xlab = "", ylab = "")
lapply(1:length(treatments),function(tr){
  Tr<-treatments[[tr]]
  axis(1, at = ACP$co[Tr,1])
  text(ACP$co[Tr,1], 0, labels = rownames(ACP$co[Tr,]),col=col[tr],cex=1.6)
})

contrib<-inertia.dudi(ACP,row.inertia=T, col.inertia=F)$row.rel
contrib_1<-contrib[order(abs(contrib[,"Axis1"]),decreasing = T),][1:15,]
contrib_2<-contrib[order(abs(contrib[,"Axis2"]),decreasing = T),][1:15,]

Sp1<-ACP$li[rownames(contrib_1),]
Sp2<-ACP$li[rownames(contrib_2),]

plot(Sp1[,"Axis1"],Sp1[,"Axis2"],axes=FALSE, type = "n", xlab = "", ylab = "")
text(Sp1[,"Axis1"], Sp1[,"Axis2"], labels = rownames(Sp1),cex=0.7)


plot(Sp2[,"Axis1"],Sp2[,"Axis2"],axes=FALSE, type = "n", xlab = "", ylab = "")
text(Sp2[,"Axis1"], Sp2[,"Axis2"], labels = rownames(Sp2),cex=0.7)

