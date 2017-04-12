# SubSubsection 2.3.2 (Function for construction of D-optimal row-column design using array exchange algorithm) 
Doptrcd.maeA<-function(trt.N,col.N,theta,nrep,itr.cvrgval)  {
  #House keeping
  arrays.1=t(combn(trt.N,2))
  arrays=rbind(arrays.1,t(rbind(arrays.1[,2],arrays.1[,1])))
  na=dim(arrays)[1]
  del.1<-matrix(10^20,na,3)
  desbest.1<-matrix(0,nrep*2,col.N)
  doptbest.1<-matrix(0,nrep,2)
  for(irep in 1:nrep){
    des<-intcrcd.mae(trt.N,col.N)
    if(trt.N==col.N&trt.N>3&irep<(trt.N-1)) {in.desns=matrix(0,(trt.N-3)*2,col.N)
    in.desns0=rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1))
    for(i in 1:(trt.N-3)) {in.desns01=cbind(rbind(seq(1,(trt.N-i)),c(seq(1,(trt.N-i))[2:(trt.N-i)],1)),rbind(rep(1,i),((trt.N-i+1):trt.N))); in.desns[c((i-1)*2+1,i*2),]=in.desns01}
    in.desns=rbind(rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1)),in.desns)
    des=in.desns[c((irep-1)*2+1,irep*2),]}
    cmat<-cmatrcd.mae(trt.N,col.N,theta,des)
    degv<-sort(eigen(cmat)$values)
    degvp<-degv[2:length(degv)]
    dopt<-prod(1/degvp) 
    dcold=dopt
    descold=t(des)
    cdel=100
    while( abs(cdel)>=0.000000001){
      i=1;
    ivaldcold={}
    for (i in 1:col.N){
      j=1;
      for (j in 1:na){
        temp=descold[i,]
        if(all(descold[i,]==arrays[j,]))  {dopt=dcold; del.1[j,]<-c(j,(dcold-dopt),dopt); next}
        descold[i,]=arrays[j,]
        trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
        R.trt<-t(trtin)%*%trtin
        if (rankMatrix(R.trt)[1]<trt.N)  {dopt=10^20; del.1[j,]<-c(j,(dcold-dopt),dopt); next}
        cmato=cmatrcd.mae(trt.N,col.N, 0,t(descold))
        egv<-sort(eigen(cmato)$values)
        if(egv[2]<0.000001) {dopt=10^20; del.1[j,]<-c(j,(dcold-dopt),dopt); next}
        cmat=cmatrcd.mae(trt.N,col.N,theta,t(descold))
        degv<-sort(eigen(cmat)$values)
        degvp<-degv[2:length(degv)]
        dopt<-prod(1/degvp) 
        del.n<-del.1[j,]<-c(j,(dcold-dopt),dopt)
        descold[i,]=temp
      }
      del.1<-del.1[order(del.1[,3]),]
      delbest=t(del.1[1,])
      descold[i,]=arrays[delbest[1],]
      dcold=delbest[3]
      cdel=delbest[2]
      ivaldcold=rbind(ivaldcold, c(i,dcold))
      if(i>itr.cvrgval) if(all(ivaldcold[c(i-(itr.cvrgval-2),i),2]==ivaldcold[i-(itr.cvrgval-1),2])) break
    }
    #print(c(000,irep,dcold,cdel,000))
  }
  #dopt0=dcold
  cdel<-1000
  while( abs(cdel)>=0.000000001){
    dopt=dcold
    #desg<-graph(t(descold))
    #plot(desg)
    del.2<-matrix(10^20,col.N+1,3)
    del.2[col.N+1,]<-c(col.N+1,0,dcold)
    for(i in 1:col.N){
      temp=descold[i,]
      descold[i,]=rev(descold[i,])
      cmato=cmatrcd.mae(trt.N,col.N, 0,t(descold))
      egv<-sort(eigen(cmato)$values)
      if(egv[2]<0.000001) {dopt2=10^20; del.2[i,]<-c(i,(dcold-dopt2),dopt2); next}
      cmat=cmatrcd.mae(trt.N,col.N,theta,t(descold))
      degv<-sort(eigen(cmat)$values)
      degvp<-degv[2:length(degv)]
      dopt2<-prod(1/degvp) 
      del.2[i,]<-c(i,(dcold-dopt2),dopt2)
      descold[i,]=temp
    }
    del.2<-del.2[order(del.2[,3]),]
    delbest=t(del.2[1,])
    if(delbest[1]<=col.N) {descold[delbest[1],]=rev(descold[delbest[1],]); cdel=delbest[2]; dcold=delbest[3]} else {cdel=0}
    #print(delbest[1]<=col.N)
    #desg<-graph(t(descold))
    #plot(desg,main=paste(dopt-dcold,sep=" / "))
    #print(del.2)
    #print(cdel)
    
    #cat("\n", dopt-dcold,"\n")
    #cdel<-dopt-dcold
    #print(cdel)
    #print(c(111,irep,dcold,dopt0-dcold,cdel,111))
  }
  #desg<-graph(t(descold))
  #plot(desg,main=paste(dopt-dcold,sep=" / "))
  #print(del.2)
  #print(cdel)
  #cat("\n", dopt-dcold,"\n")
  #"============================================================="
  
  if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    doptbest.1[irep,]=c(irep,dcold)
  }
  best=doptbest.1[order(doptbest.1[,2]),]
  nb=best[1,1]
  Dscore<-best[1,2]
  Doptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))
  cnames=paste0("Ary",1:col.N)
  dimnames(Doptde)=list(c("Dye 1:", "Dye 2:"),cnames)
  Dopt_sum2<-list("v"=trt.N,"b"=col.N,theta=theta,nrep=nrep,itr.cvrgval=itr.cvrgval, "OptdesF"=Doptde,"Optcrtsv" =Dscore)
  return(Dopt_sum2)
}#End of SubSubsection 2.3.2 (Doptrcd.maeA function) construction of D-optimal row-column design using array exchange algorithm
#End of Subsection 2.3 (Function for construction of D-optimal row-column design)    
