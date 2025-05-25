# SubSubsection 2.2.2 (Function for construction of MV-optimal row-column design using array exchange algorithm) 
MVoptrcd.maeA<-function(trt.N,col.N,theta,nrep,itr.cvrgval) {
  #House keeping
  arrays.1=t(combn(trt.N,2))
  arrays=rbind(arrays.1,t(rbind(arrays.1[,2],arrays.1[,1])))
  na=dim(arrays)[1]
  ii=2
  trco=cbind(matrix(1,trt.N-1),-diag(1,trt.N-1,trt.N-1))
  while(ii<=trt.N-1){
    if (ii==trt.N-1){
      trco1=cbind(matrix(0,1,trt.N-2),matrix(1,trt.N-ii),-diag(1,trt.N-ii,trt.N-ii))}
    else
    {trco1=cbind(matrix(0,trt.N-ii,trt.N-(trt.N-ii+1)),matrix(1,trt.N-ii),-diag(1,trt.N-ii,trt.N-ii))}
    trco=rbind(trco,trco1)
    ii=ii+1
  }
  del.1<-matrix(10^20,na,3)
  desbest.1<-matrix(0,nrep*2,col.N)
  MVoptbest.1<-matrix(0,nrep,2)
  for(irep in 1:nrep){
    des<-intcrcd.mae(trt.N, col.N)
    if(trt.N==col.N&trt.N>3&irep<(trt.N-1)) {in.desns=matrix(0,(trt.N-3)*2,col.N)
    in.desns0=rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1))
    for(i in 1:(trt.N-3)) {in.desns01=cbind(rbind(seq(1,(trt.N-i)),c(seq(1,(trt.N-i))[2:(trt.N-i)],1)), rbind(rep(1,i),((trt.N-i+1):trt.N))); in.desns[c((i-1)*2+1,i*2),]=in.desns01}
    in.desns=rbind(rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1)),in.desns)
    des=in.desns[c((irep-1)*2+1,irep*2),]}
    cmat<-cmatrcd.mae(trt.N,col.N,theta,des)
    invc=ginv(cmat)
    invcp=trco%*%invc%*%t(trco);
    MVopt =max(diag(invcp));
    MVcold=MVopt
    descold=t(des)
    cdel=100
    while( abs(cdel)>=0.000000001){
      ivalMVcold={}
    for (i in 1:col.N){
      j=1;
      for (j in 1:na){
        temp=descold[i,]
        if(all(descold[i,]==arrays[j,]))  {MVopt=MVcold; del.1[j,]<-c(j,(MVcold-MVopt),MVopt); next}
        descold[i,]=arrays[j,]
        trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
        R.trt<-t(trtin)%*%trtin
        if (rankMatrix(R.trt)[1]<trt.N)  {MVopt=10^20; del.1[j,]<-c(j,(MVcold-MVopt),MVopt); next}
        cmato=cmatrcd.mae(trt.N,col.N, 0,t(descold))
        egv<-sort(eigen(cmato)$values)
        if(egv[2]<0.000001) {MVopt=10^20; del.1[j,]<-c(j,(MVcold-MVopt),MVopt); next}
        cmat=cmatrcd.mae(trt.N,col.N,theta,t(descold))
        invc=ginv(cmat)
        invcp=trco%*%invc%*%t(trco);
        MVopt =max(diag(invcp));
        del.n<-del.1[j,]<-c(j,(MVcold-MVopt),MVopt)
        descold[i,]=temp
      }
      del.1<-del.1[order(del.1[,3]),]
      delbest=t(del.1[1,])
      descold[i,]=arrays[delbest[1],]
      MVcold=delbest[3]
      cdel=delbest[2]
      ivalMVcold=rbind(ivalMVcold, c(i,MVcold))
      if(i>itr.cvrgval) if(all(ivalMVcold[c(i-(itr.cvrgval-2),i),2]==ivalMVcold[i-(itr.cvrgval-1),2])) break
    }
    #print(c(000,irep,MVcold,cdel,000))
  }
  #MVopt0=MVcold
  cdel<-1000
  while( abs(cdel)>=0.000000001){
    MVopt=MVcold
    #desg<-graph(t(descold))
    #plot(desg)
    del.2<-matrix(10^20,col.N+1,3)
    del.2[col.N+1,]<-c(col.N+1,0,MVcold)
    for(i in 1:col.N){
      temp=descold[i,]
      descold[i,]=rev(descold[i,])
      cmato=cmatrcd.mae(trt.N,col.N, 0,t(descold))
      egv<-sort(eigen(cmato)$values)
      if(egv[2]<0.000001) {MVopt2=10^20; del.2[i,]<-c(i,(MVcold-MVopt2),MVopt2); next}
      cmat=cmatrcd.mae(trt.N,col.N,theta,t(descold))
      MVopt2=sum(diag(ginv(cmat)))
      del.2[i,]<-c(i,(MVcold-MVopt2),MVopt2)
      descold[i,]=temp
    }
    del.2<-del.2[order(del.2[,3]),]
    delbest=t(del.2[1,])
    if(delbest[1]<=col.N) {descold[delbest[1],]=rev(descold[delbest[1],]); cdel=delbest[2]; MVcold=delbest[3]} else {cdel=0}
    #print(delbest[1]<=col.N)
    #desg<-graph(t(descold))
    #plot(desg,main=paste(MVopt-MVcold,sep=" / "))
    #print(del.2)
    #print(cdel)
    
    #cat("\n", MVopt-MVcold,"\n")
    #cdel<-MVopt-MVcold
    #print(cdel)
    #print(c(111,irep,MVcold,MVopt0-MVcold,cdel,111))
  }
  #desg<-graph(t(descold))
  #plot(desg,main=paste(MVopt-MVcold,sep=" / "))
  #print(del.2)
  #print(cdel)
  #cat("\n", MVopt-MVcold,"\n")
  #"============================================================="
  
  if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    MVoptbest.1[irep,]=c(irep,MVcold)
  }
  best=MVoptbest.1[order(MVoptbest.1[,2]),]
  nb=best[1,1]
  MVscore<-best[1,2]
  MVoptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  if(trt.N!=3) {tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))}
  cnames=paste0("Ary",1:col.N)
  dimnames(MVoptde)=list(c("Dye 1:", "Dye 2:"),cnames)
  MVopt_sum2<-list("v"=trt.N,"b"=col.N,theta=theta,nrep=nrep,itr.cvrgval=itr.cvrgval, "OptdesF"=MVoptde,"Optcrtsv" =MVscore)
  return(MVopt_sum2)
}#End of SubSubsection 2.2.2 (MVoptrcd.maeA function) construction of MV-optimal row-column design using array exchange algorithm
#End of Subsection 2.2 (Function for construction of MV-optimal row-column design)    

