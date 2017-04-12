\name{graphoptrcd.mae}
\alias{graphoptrcd.mae}
\title{
Creates the graphical layout of resultant optimal row-column design
}
\description{
Creates the graphical layout of resultant A-, MV-, D- or E-optimal or near-optimal row-column design on a separate pop-up GUI tcltk window.
}
\usage{
graphoptrcd.mae(trt.N, col.N, theta, OptdesF, Optcrit, cbVal2)
}
\arguments{
  \item{trt.N}{
integer, specifying  number of treatments, \code{v}. 
}
  \item{col.N}{
integer, specifying  number of arrays, \code{b}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{OptdesF}{
matrix, a \code{2 x b} obtained optimal or near-optimal row-column design.
}
  \item{Optcrit}{
character, specifying the optimality criteria to be used. Thus, \code{Optcrit} takes the letter \code{"A"}, \code{"MV"}, \code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal row-column designs, respectively.
}
  \item{cbVal2}{
checkbox value.  It takes a value of zero or one. The default value of \code{cbVal2} is 0. Thus, if \code{cbVal2 = 0}, the function will display the graphical layout of the resultant optimal design generated using the treatment exchange algorithm. Similarly, if \code{cbVal2 = 1}, the function will display the graphical layout of the resultant optimal design generated using the array exchange algorithm. 
}
}
\details{
Detail discussions  concerning the constructions of a graphs can be found in \code{\link{igraph}} R package.
}
\value{
Returns the graphical layout of the resultant optimal or near-optimal row-column design '\code{OptdesF}' on a separate pop-up window. Furthermore, the function  \code{graphoptrcd.mae} saves the graphical layout of the resultant optimal or near-optimal row-column design in .pdf format in a working subdirectory. 
}
\author{
Legesse Kassa Debusho, Dibaba Bayisa Gemechu, and Linda Haines
}
\note{
When closing a pop-up window for graphical layout of the resultant designs (Graph plot),  if  the window is closed by 
clicking on the red button with "X" sign (top-right), the warning message  "Warning message: In \code{rm(list = cmd,} \code{envir = .tkplot.env):} object 'tkp \dots' not found"
will occur in R console irrespective of what command is executed next. To resolve this warning message, click 
on "\code{close}" menu that is located at the top-left of the graph plot pop-up window when closing this window. 
 }

\seealso{
\code{\link{igraph}}
}
\examples{
    ##To create the graphical layout of the D-optimal row-column design
    ##obtained using the treatment exchange algorithm for 
    
    trt.N <- 10  #Number of treatments
    
    col.N <- 10  #Number of arrays
    
    theta <- 0.2 #theta value 
    
    OptdesF <- rbind(1:10, c(2:10, 1)) #D-optimal design (loop design)

    graphoptrcd.mae(trt.N = 10, col.N = 10, theta = 0.2, OptdesF, Optcrit = "D", cbVal2 = 0)
}
\keyword{Graphical layout}
