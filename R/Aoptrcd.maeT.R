#Subsection 2.1: Function for search of A-optimal or near-optimal row-column designs
# SubSubsection 2.1.1 (Function for construction of A-optimal row-column designs using treatment exchange algorithm) 
Aoptrcd.maeT<-function(trt.N,col.N,theta,nrep,itr.cvrgval) {
  #House keeping
  del.1<-matrix(1000,trt.N,3)
  desbest.1<-matrix(0,nrep*2,col.N)
  aoptbest.1<-matrix(0,nrep,2)
  #Start iteration
  for(irep in 1:nrep){
    #Initial design with its corresponding Ascore value
    des<-intcrcd.mae(trt.N,col.N)
    if(trt.N==col.N&trt.N>3&irep<(trt.N-1)) {in.desns=matrix(0,(trt.N-3)*2,col.N)
    in.desns0=rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1))
    for(i in 1:(trt.N-3)) {in.desns01=cbind(rbind(seq(1,(trt.N-i)),c(seq(1,(trt.N-i))[2:(trt.N-i)],1)), rbind(rep(1,i),((trt.N-i+1):trt.N))); in.desns[c((i-1)*2+1,i*2),]=in.desns01}
    in.desns=rbind(rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1)),in.desns)
    des=in.desns[c((irep-1)*2+1,irep*2),]}
    cmat<-cmatrcd.mae(trt.N,col.N,theta,des)
    aopt=sum(diag(ginv(cmat)))
    acold=aopt
    descold=t(des)
    #deletion difference
    cdel=100
    while( abs(cdel)>=0.000000001){
    i=1;
    ivalacold={}
    for(i in 1:col.N){
      for (m in 1:2){
        j=1;
        for(j in 1:trt.N){
          temp=descold[i,]
          if(m==1) {
            if(j==descold[i,1]|j==descold[i,2]) {aopt=acold; del.1[j,]<-c(descold[i,1],(acold-aopt),aopt); next} else { descold[i,]=c(j,descold[i,2])}}
          if(m==2) {
            if(descold[i,2]==j|j==descold[i,1]) {aopt=acold; del.1[j,]<-c(descold[i,2],(acold-aopt),aopt); next} else { descold[i,]=c(descold[i,1],j)}}
          
          trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
          R.trt<-t(trtin)%*%trtin
          if (rankMatrix(R.trt)[1]<trt.N)  {aopt=acold; descold[i,]=temp; if(m==1) {del.1[j,]<-c(descold[i,1],(acold-aopt),aopt)} else {
           del.1[j,]<-c(descold[i,2],(acold-aopt),aopt)}; next}
          cmato=cmatrcd.mae(trt.N,col.N, 0,t(descold))
          egv<-sort(eigen(cmato)$values)
          if(egv[2]<0.000001) {aopt=acold; descold[i,]=temp; if(m==1){del.1[j,]<-c(descold[i,1],(acold-aopt),aopt)} else {
           del.1[j,]<-c(descold[i,2],(acold-aopt),aopt)}; next}
          cmat=cmatrcd.mae(trt.N,col.N,theta,t(descold))
          aopt=sum(diag(ginv(cmat)))
          del.n<-del.1[j,]<-c(j,(acold-aopt),aopt)
          descold[i,]=temp
        }
        del.1<-del.1[order(del.1[,3]),]
        delbest=t(del.1[1,])
        if (m==1) {
          if (delbest[1]==descold[i,2]) {descold[i,]=descold[i,]}  else 
          {descold[i,]=c(delbest[1],descold[i,2]); cdel=delbest[2]; acold=delbest[3]}} else {
            if (descold[i,1]==delbest[1]) {descold[i,]= descold[i,]} else 
            {descold[i,]=c(descold[i,1],delbest[1]); cdel=delbest[2]; acold=delbest[3]}}
      }
      ivalacold=rbind(ivalacold, c(i,acold))
      if(i>itr.cvrgval) if(all(ivalacold[c(i-(itr.cvrgval-2),i),2]==ivalacold[i-(itr.cvrgval-1),2])) break
    }
    #print(c(000,irep,acold,cdel,000))
}
    #aopt0=acold
    cdel<-1000
    while( abs(cdel)>=0.000000001){
      aopt=acold
      #desg<-graph(t(descold))
      #plot(desg)
      del.2<-matrix(1000,col.N+1,3)
      del.2[col.N+1,]<-c(col.N+1,0,acold)
      for(i in 1:col.N){
        temp=descold[i,]
        descold[i,]=rev(descold[i,])
        cmato=cmatrcd.mae(trt.N,col.N, 0,t(descold))
        egv<-sort(eigen(cmato)$values)
        if(egv[2]<0.000001) {aopt2=1000; del.2[i,]<-c(i,(acold-aopt2),aopt2); next}
        cmat=cmatrcd.mae(trt.N,col.N,theta,t(descold))
        aopt2=sum(diag(ginv(cmat)))
        del.2[i,]<-c(i,(acold-aopt2),aopt2)
        descold[i,]=temp
      }
      del.2<-del.2[order(del.2[,3]),]
      delbest=t(del.2[1,])
      if(delbest[1]<=col.N) {descold[delbest[1],]=rev(descold[delbest[1],]); cdel=delbest[2]; acold=delbest[3]} else {cdel=0}
      #print(delbest[1]<=col.N)
      #desg<-graph(t(descold))
      #plot(desg,main=paste(aopt-acold,sep=" / "))
      #print(del.2)
      #print(cdel)
      
      #cat("\n", aopt-acold,"\n")
      #cdel<-aopt-acold
      ##print(cdel)
      #print(c(111,irep,acold,aopt0-acold,cdel,111))
    }
    #desg<-graph(t(descold))
    #plot(desg,main=paste(aopt-acold,sep=" / "))
    #print(del.2)
    #print(cdel)
    #cat("\n", aopt-acold,"\n")
    #"============================================================="
    
    next.it<- if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    aoptbest.1[irep,]=c(irep,acold)
    #print(c(222,irep,acold,aopt0-acold,222))
  }
  
  best=aoptbest.1[order(aoptbest.1[,2]),]
  #print(best)
  nb=best[1,1]
  Ascore<-best[1,2]
  Aoptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  if(trt.N!=3) {tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))}
  cnames=paste0("Ary",1:col.N)
  dimnames(Aoptde)=list(c("Dye 1:", "Dye 2:"),cnames)
  Aopt_sum2<-list("v"=trt.N,"b"=col.N,theta=theta,nrep=nrep,itr.cvrgval=itr.cvrgval, "OptdesF"=Aoptde,"Optcrtsv" =Ascore)
  return(Aopt_sum2)
}#End of SubSubsection 2.1.1 (Aoptrcd.maeT function) construction of A-optimal row-column design using treatment exchange algorithm
