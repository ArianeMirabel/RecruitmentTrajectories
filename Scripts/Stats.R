load("DB/Turnover_toInit")
treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

Maxs<-cbind(apply(Turn,1,max),treats[,"treat"])
colnames(Maxs)<-c("Max","treat")

cor(Maxs[,"Max"],Maxs[,"treat"],method="spearman")
