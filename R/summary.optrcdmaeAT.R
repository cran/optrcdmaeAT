summary.optrcdmaeAT<-function(object,...)
{
optrcd_maes<-object
optrcd_maes$grphlt<-make_graph(as.numeric(as.factor(object$OptdesF)), directed = TRUE)
class(optrcd_maes)<-"summary.optrcdmaeAT"
optrcd_maes
}
