#Section 3 computation of the information matrix (C-matrix) for a given row-column design (des) 
cmatrcd.mae<-function(trt.N,col.N,theta,des){
  k=2
  trtin<-contrasts(as.factor(des),contrasts=FALSE)[as.factor(des),]
  col.1<-rep(1:col.N,each=2)
  colin<-contrasts(as.factor(col.1),contrasts=FALSE)[as.factor(col.1),]
  vec.1<-rep(1,col.N*2)
  R.trt<-t(trtin)%*%trtin
  N.tb<-t(trtin)%*%colin
  r.trt<-t(trtin)%*%vec.1
  forrow.1=contrasts(as.factor(c(des[1,],seq(1:trt.N))),contrasts=FALSE)[as.factor(c(des[1,],seq(1:trt.N))),]
  forrow.2=contrasts(as.factor(c(des[2,],seq(1:trt.N))),contrasts=FALSE)[as.factor(c(des[2,],seq(1:trt.N))),]
  mmmat.r=t(rbind(t(as.matrix(colSums(forrow.1)-1)), t(as.matrix(colSums(forrow.2)-1))))%*% 
    rbind(t(as.matrix(colSums(forrow.1)-1)), t(as.matrix(colSums(forrow.2)-1)))
  cmat<-R.trt-(1/k)*(N.tb%*%t(N.tb))-(1/col.N)*(mmmat.r)+(1/(col.N*k))*(r.trt%*%t(r.trt))+
    theta*((1/k)*(N.tb%*%t(N.tb))-(1/(col.N*k))*(r.trt%*%t(r.trt)))
  #list(cmat0=cmatc,cmat=cmat)
  cmat
}#End of Section 3 (computation of C-matrix, cmatrcd.mae)
