print.summary.optrcdmaeAT<-function(x,...)
{
cat("\n        ---------------------------------------    \n")
cat("Title:  ",x$Optcrit,"-optimal or near-optimal row-column design        ","Date: ", format(Sys.time(), "%a %b %d %Y %H:%M:%S"),"\n",sep="")
cat("        ---------------------------------------    \n")
cat("Call:\n")
print(x$call)
cat("\nparametric combinations:\n")
cat("\nNumber of treatments:          ", x$v, "\n")
cat("Number of arrays:              ", x$b, "\n")
cat("Theta value:                   ",  x$theta, "\n")
cat("Number of replications:        ",  x$nrep, "\n")
cat("Number of exchange iteration:  ", x$itr.cvrgval, "\n")
cat("Algorithm used:                ", x$Alg, "\n")
cat("Optimality criterion used:      ", x$Optcrit, "-optimality criteria\n",sep="")
cat("\nResultant ",x$Optcrit,"-optimal or near-optimal row-column design:\n",sep="")
cat("\n")
print(data.frame(x$OptdesF))
cat("\n")
cat(x$Optcrit,"-Score value:  ", x$Optcrtsv, "\n",sep="")
plot(x$grphlt,edge.arrow.size=0.5, vertex.size=28, margin=0.5,
       layout=layout.kamada.kawai,vertex.color="cyan",edge.color="black")
title(paste("Graphical layout of", paste(x$Optcrit,"-optimal or near-optimal row-column design",sep=""),sep=" "), 
      sub = NULL,cex.main = 1,   font.main= 1, col.main= "black")
mtext(paste("using"," ",x$Alg," algorithm for:",sep=""), line = 0.5, col = "black", font = 1)
mtext(paste("(v, b, theta) =", " (",paste(x$v, x$b, x$theta, sep=", "),")",sep=""), line = -0.50, col = "blue", font = 1)
cat("\n", x$file_loc2,"\n", x$file_loc,"\n")
cat("\n")
}
