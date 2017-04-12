# Section 1
# General functions & TCL/TK functions
fixparrcd.mae<-function(Optcrit){
  trt.I<-tclVar(paste("3"));#Number of treatments (default)
  blk.I<-tclVar(paste("3"));#Number of arrays (default)
  theta.I<-tclVar(paste("0.0"));#Theta value (default)
  rep.I<-tclVar(paste("10"));#Number of replications (default)
  itrcvrgval.I<-tclVar(paste("6"));#Initial convergence value (default)
  cbValue.I<-tclVar("0")
  cbValue2.I<-tclVar("0")  
  cbValue3.I<-tclVar("0")  
#"=============================================================================="
optcrtF<-function(Optcrit){
   nrep<-as.numeric(tclvalue(rep.I))
   trt.N<-as.numeric(tclvalue(trt.I))
   col.N<-as.numeric(tclvalue(blk.I))
   theta<-as.numeric(tclvalue(theta.I))
   cbVal<-as.numeric(tclvalue(cbValue.I))
   itr.cvrgval<-as.numeric(tclvalue(itrcvrgval.I))
   if(itr.cvrgval>col.N) itr.cvrgval<-col.N
  cbVal2<-as.numeric(tclvalue(cbValue2.I))
  cbVal3<-as.numeric(tclvalue(cbValue3.I))
if(Optcrit=="A") if(cbVal2==0) {optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="A",Alg="trtE")} else {
   optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="A",Alg="arrayE")} 
if(Optcrit=="MV") if(cbVal2==0) {optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="MV",Alg="trtE")} else {
   optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="MV",Alg="arrayE")} 
if(Optcrit=="D") if(cbVal2==0) {optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="D",Alg="trtE")} else {
   optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="D",Alg="arrayE")}
if(Optcrit=="E") if(cbVal2==0) {optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="E",Alg="trtE")} else {
   optrcdFS=optrcdmaeAT(trt.N,col.N,theta,nrep,itr.cvrgval,Optcrit="E",Alg="arrayE")}
if(cbVal3=="1") optrcdFS=summary(optrcdFS) 
print(optrcdFS)
if(cbVal=="1") {graphoptrcd.mae(trt.N, col.N,theta,optrcdFS$OptdesF,Optcrit,cbVal2)}
}
#"=============================================================================="

  fiPar<-tktoplevel()
  fontHeading<- tkfont.create(family="times",size=40,weight="bold")
  fontHeading3<-tkfont.create(family="times",size=10,weight="bold")
  fontHeading2<-tkfont.create(family="times",size=14,weight="bold")
  fontHeading4<-tkfont.create(family="times",size=12,weight="bold")
  tkwm.title(fiPar,"Set parameter values")
  fiParF<-tkframe(fiPar)
  fiParFup<- tkframe(fiParF,relief="groove",borderwidth=2)
  fiParFmid<-tkframe(fiParF,relief="sunken",borderwidth=2)
  fiParFlow<-tkframe(fiParF,relief="sunken",borderwidth=2)
  fiParFlow2<-tkframe(fiParF,relief="groove",borderwidth=2)
  trt.in<-trt.I
  blk.in<-blk.I
  theta.in<-theta.I
  rep.in<-rep.I
  itrcvrgval.in<-itrcvrgval.I
  cbValue.in<-cbValue.I
  cbValue2.in<-cbValue2.I
  cbValue3.in<-cbValue3.I
  trt.i1<-tkentry(fiParFup,width=11,textvariable=trt.in)
  blk.i1<-tkentry(fiParFup,width=11,textvariable=blk.in)
  theta.i1<-tkentry(fiParFup,width=11,textvariable=theta.in)
  rep.i1<-tkentry(fiParFup,width=11,textvariable=rep.in)
  itrcvrgval.i1<-tkentry(fiParFup,width=11,textvariable=itrcvrgval.in)
  cbValue.i1<-tkcheckbutton(fiPar)
  tkconfigure(cbValue.i1,variable=cbValue.in)
  cbValue.i3<-tkcheckbutton(fiPar)
  tkconfigure(cbValue.i3,variable=cbValue3.in)
  TrtEx<-tkradiobutton(fiPar)
  ArEx<-tkradiobutton(fiPar)
  tkconfigure(TrtEx,variable=cbValue2.in,value="0")
  tkconfigure(ArEx,variable=cbValue2.in,value="1")
  tkgrid(tklabel(fiParF,text="       Fix Value of Parameters     ",font=fontHeading2))
  tkgrid(tklabel(fiParFup,text="Number of treatments:                   "),trt.i1)
  tkgrid(tklabel(fiParFup,text="Number of arrays:                           "),blk.i1)
  tkgrid(tklabel(fiParFup,text="Theta value:                                       "),theta.i1)
  tkgrid(tklabel(fiParFup,text="Number of replications:                  "),rep.i1)
  tkgrid(tklabel(fiParFup,text="Iterations for exchange Procedure: "),itrcvrgval.i1)
  tkgrid(tklabel(fiParFlow,text="Algorithm to be used:      ", font=fontHeading4))
  tkgrid(tklabel(fiParFlow,text="   Treatment Exchange  "),TrtEx)
  tkgrid(tklabel(fiParFlow,text="Array Exchange       "),ArEx)
  tkgrid(tklabel(fiParFmid,text="Show graphical layout        ", font=fontHeading4),cbValue.i1)
  tkgrid(tklabel(fiParFmid,text="Show Summary result        ", font=fontHeading4),cbValue.i3)
  tkgrid(tklabel(fiParFlow,text=" ",font=0.01))#empty line
  exitFP<-function() {
    closeQ=tkmessageBox(title = "Exit set parameter values", message = "You are leaving set parameter values window",
                        icon = "info", type = "okcancel", default = "cancel")
    if(as.character(closeQ)=="ok") tkdestroy(fiPar)
  }#end of exitFP
  exit.but2a<-tkbutton(fiParF,text="Exit",command=exitFP,width=9)
  serch.but<-function(Optcrit){
    if (Optcrit=="A") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search A-opt row-column design",command=function()optcrtF("A"),width=28), exit.but2a)
    if (Optcrit=="MV") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search MV-opt row-column design",command=function()optcrtF("MV"),width=29), exit.but2a)
    if (Optcrit=="D") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search D-opt row-column design",command=function()optcrtF("D"),width=28), exit.but2a)
    if (Optcrit=="E") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search E-opt row-column design",command=function()optcrtF("E"),width=28), exit.but2a)
    return(but.A)
  }#end of serch.but
  tkgrid(tklabel(fiParFlow2,text="                                                    ",font=0.00001))
  serch.but(Optcrit)
  tkgrid(tklabel(fiParFlow2,text="                                             ",font=0.00001))
  tkgrid(fiParFup)
  tkgrid(tklabel(fiParF,text="---------------------------------", font=fontHeading3))
  tkgrid(fiParFlow)
  tkgrid(tklabel(fiParF,text="---------------------------------", font=fontHeading3))
  tkgrid(fiParFmid)
  tkgrid(tklabel(fiParF,text="---------------------------------", font=fontHeading3))
  tkgrid(fiParFlow2)
  tkgrid(fiParF) 
}#end of fixparrcd.mae
# End of Section 1 (general functions & TCL/TK functions)

