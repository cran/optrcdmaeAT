optrcdmaeAT.default<-function(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="",Alg="",...)
{
  trt.N=as.numeric(trt.N)
  col.N=as.numeric(col.N)
  theta=as.numeric(theta)
  nrep=as.numeric(nrep)
  itr.cvrgval=as.numeric(itr.cvrgval)
  #"===================================================================================================="
  if(is.na(theta)|theta<0|theta>1){
    tkmessageBox(title="Error",message=paste("Please insert correct value of theta, it should be between 0 and 1 inclusive of 0 and 1, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'theta', it should be between 0 and 1, inclusive of 0 and 1")
  }#end of if
  if(is.na(trt.N)|is.na(col.N)|trt.N!=round(trt.N)|col.N!=round(col.N)) {
    tkmessageBox(title="Error",message=paste("Please insert correct format of the number of treatments and arrays. The number of treatments and arrays should be an integer, click OK to reset the values.",sep=""),icon = "error"); 
    stop("Wrong format of 'trt.N' and/or 'col.N', both should be an integer")
  }#end of if
  if(trt.N<3|col.N<3){ 
    tkmessageBox(title="Error",message=paste("The number of arrays and treatments should be greater than or equal to 3, click Ok to reset.",sep=""),icon = "error"); 
    stop("Very small value of number of treatments and/or arrays, minimum value of the two is 3")
  }#end of if
  if(trt.N-col.N>0){ 
    tkmessageBox(title="Error",message=paste("The number of arrays should be greater than or equal to the number of treatments, click Ok to reset.",sep=""),icon = "error"); 
    stop("The number of treatments are larger than the number of arrays, 'trt.N' should be less than or equal to 'col.N' ")
  }#end of if
  if(trt.N>10|col.N>10){ 
    tkmessageBox(title="Information",message=paste("This might take some minutes, please be patient...",sep=""))
  }#end of if
  if(is.na(itr.cvrgval)|itr.cvrgval<2|itr.cvrgval!=round(itr.cvrgval)){
    tkmessageBox(title="Error",message=paste("The number of iteration for the exchange procedure should be a positive integer greater than or equal to two, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'itr.cvrgval', it should be greater than or equal to two (only positive integer values)")
  }#end of if
  if(is.na(nrep)|nrep<2|nrep!=round(nrep)){
    tkmessageBox(title="Error",message=paste("The number of replications should be a positive integer greater than or equal to two, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'nrep', it should be greater than or equal to two (only positive integer values)")
  }#end of if
  #"===================================================================================================="
  if(!is.element(Optcrit,c("A","MV","D","E"))){stop("The optimality criterion 'Optcrit' is not correctly specified")}
  if(!is.element(Alg,c("trtE","arrayE"))){stop("The algorithm 'Alg' is not correctly specified")}
  if(itr.cvrgval>col.N) itr.cvrgval<-col.N
  if(Alg=="trtE") {
    if(Optcrit=="A") {
      optrcd_mae<-Aoptrcd.maeT(trt.N,col.N,theta,nrep,itr.cvrgval)} else if(
        Optcrit=="MV") {
        optrcd_mae<-MVoptrcd.maeT(trt.N,col.N,theta,nrep,itr.cvrgval)} else if(
          Optcrit=="D") {
          optrcd_mae<-Doptrcd.maeT(trt.N,col.N,theta,nrep,itr.cvrgval)} else if(
            Optcrit=="E") {
            optrcd_mae<-Eoptrcd.maeT(trt.N,col.N,theta,nrep,itr.cvrgval)} else{
              stop("The optimality criterion is not specified")}
    optrcd_mae$Alg="Treatment exchange"} else if(Alg=="arrayE") {
      if(Optcrit=="A") {
        optrcd_mae<-Aoptrcd.maeA(trt.N,col.N,theta,nrep,itr.cvrgval)} else if(
          Optcrit=="MV") {
          optrcd_mae<-MVoptrcd.maeA(trt.N,col.N,theta,nrep,itr.cvrgval)} else if(
            Optcrit=="D") {
            optrcd_mae<-Doptrcd.maeA(trt.N,col.N,theta,nrep,itr.cvrgval)} else if(
              Optcrit=="E") {
              optrcd_mae<-Eoptrcd.maeA(trt.N,col.N,theta,nrep,itr.cvrgval)} else{
                stop("The optimality criterion is not specified")}
      optrcd_mae$Alg="Array exchange"
    } else {stop("The algorithm is not specified")}#end of if
  optrcd_mae$call<-match.call()
  optrcd_mae$Optcrit<-Optcrit
  optrcd_mae$Cmat<-cmatrcd.mae(optrcd_mae$v,optrcd_mae$b,optrcd_mae$theta,optrcd_mae$OptdesF)
  trtin <- contrasts(as.factor(optrcd_mae$OptdesF), contrasts = FALSE)[as.factor(optrcd_mae$OptdesF), ]
  vec1 <- rep(1, optrcd_mae$b * 2)
  vec_trtr <- t(trtin) %*% vec1
  optrcd_mae$equireplicate<-all(vec_trtr==vec_trtr[1])
  optrcd_mae$vtrtrep<-t(vec_trtr)
  
  #"======================================================================================"
  titleoptrcd<-list(c("      --------------------------------------- ",paste("Title: ",Optcrit,"-optimal or near-optimal row-column design          Date:", format(Sys.time(), "%a %b %d %Y %H:%M:%S"),sep=""),
                      "      --------------------------------------- "))
  write.table(titleoptrcd, file = file.path(tempdir(), paste(Optcrit,"optrcd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  parcomb<-list(c("     Parametric combination:", "Number of treatments:", "Number of arrays:", 
                  "Theta value:", "Number of replications:","Number of exchange iteration:","Algorithm used:", "OPtimality criterion used:"," ","Design obtained:"),
                c(" ",optrcd_mae$v,optrcd_mae$b,optrcd_mae$theta,optrcd_mae$nrep,optrcd_mae$itr.cvrgval,optrcd_mae$Alg,paste(Optcrit,"-optimality criterion",sep="")," "," "))
  write.table(parcomb, file =  file.path(tempdir(), paste(Optcrit,"optrcd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  
  optde<-list("",cbind(c(" ", "Dye 1:", "Dye 2:"),rbind(paste0("Ary",1:optrcd_mae$b),optrcd_mae$OptdesF)))
  write.table(optde, file = file.path(tempdir(), paste(Optcrit,"optrcd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  write.table(list(c("",paste(Optcrit,"-score value:",sep=""), "Equireplicate:",""),c("",optrcd_mae$Optcrtsv,optrcd_mae$equireplicate,"")), file =   file.path(tempdir(), paste(Optcrit,"optrcd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  write.table(list(c("Treatment:", "Treatment replication:"),rbind(1:optrcd_mae$v,optrcd_mae$vtrtrep)), file = file.path(tempdir(), paste(Optcrit,"optrcd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
 
  optrcd_mae$file_loc<-file.path(tempdir(),  paste(Optcrit,"optrcd_",Alg,"_summary.csv",sep = ""))
  optrcd_mae$file_loc2<-paste("Summary of obtained ",Optcrit,"-optimal or near-optimal row-column design is also saved at:",sep="")
  #"======================================================================================"
  class(optrcd_mae)<-"optrcdmaeAT"
  optrcd_mae
}