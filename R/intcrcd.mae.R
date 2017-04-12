#Section 4: Function for generation of initial row-column design for given number of treatments (trt.N) and arrays (col.N)
intcrcd.mae<-function(trt.N,col.N){
  arrays.1=t(combn(trt.N,2))
  arrays=rbind(arrays.1,t(rbind(arrays.1[,2],arrays.1[,1])))
  na=dim(arrays)[1]
  Con.egv.check=0.00000001
  while(Con.egv.check<0.000001){
    All.trt.check=trt.N-1
    while(All.trt.check <trt.N){
      des.p<-if(col.N>na) {c(sample(1:na,na,replace=FALSE), sample(1:na,col.N-na,replace=TRUE))} else {
        sample(1:na,col.N,replace=FALSE)}
      des<-t(arrays[des.p,])
      trtin<-contrasts(as.factor(des),contrasts=FALSE)[as.factor(des),]
      R.trt<-t(trtin)%*%trtin
      All.trt.check<-rankMatrix(R.trt)
    }
    cmato=cmatrcd.mae(trt.N,col.N, 0,des)
    egv<-sort(eigen(cmato)$values)
    Con.egv.check<-egv[2]
  }
  return(des)
}#End of Section 4 (generation of initial design, intcrcd.mae)

