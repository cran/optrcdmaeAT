\name{optrcdmaeAT}
\alias{optrcdmaeAT}
\alias{optrcdmaeAT.default}
\alias{print.optrcdmaeAT}
\alias{summary.optrcdmaeAT}
\alias{print.summary.optrcdmaeAT}
\title{
Optimal row-column designs for two-colour cDNA microarray experiments
}
\description{
Used to compute A-, MV-, D- or E-optimal or near-optimal row-column designs for two-colour cDNA microarray experiments under either the linear fixed effects model or the linear mixed effects model settings using either the array exchange or treatment exchange algorithms of Debusho, Gemechu and Haines (2016) after adjusting to the row-column setup.}
\usage{
optrcdmaeAT(trt.N, col.N, theta, nrep, itr.cvrgval, Optcrit = "", Alg = "", ...)

\method{optrcdmaeAT}{default}(trt.N, col.N, theta, nrep, itr.cvrgval, Optcrit = "", Alg = "", ...)
\method{print}{optrcdmaeAT}(x, ...)
\method{summary}{optrcdmaeAT}(object, ...)
}
\arguments{
  \item{trt.N}{
integer, specifying number of treatments, \code{v}. 
}
  \item{col.N}{
integer, specifying number of arrays, \code{b}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{nrep}{
integer, specifying number of replications of the optimization procedure. 
}
  \item{itr.cvrgval}{
integer, specifying number of iterations required for convergence during the exchange procedure. 
}
  \item{Optcrit}{
character, specifying the optimality criteria to be used. \code{Optcrit} takes the letter \code{"A"}, \code{"MV"}, \code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal row-column designs, respectively.
}
  \item{x}{
the object to be printed.
}
  \item{object}{
an object of class \code{"optrcdmaeAT"}.
}
  \item{Alg}{
character string used to specify the algorithm to be used. Possible values of \code{Alg} are \code{Alg="trtE"} for the treatment exchange algorithm and \code{Alg="arrayE"} for the array exchange algorithm: see 'Details'.
}
  \item{\dots}{
not used.
}
}
\details{
\code{optrcdmaeAT} computes optimal or near-optimal row-column design for the two-colour cDNA microarray experiments  
where the interest is in a comparison of all possible elementary treatment contrasts. The function computes A-, MV-, D- and E-optimal 
or near optimal row-column designs via calling of eight sub-functions \code{\link{Aoptrcd.maeT}}, \code{\link{Aoptrcd.maeA}}, 
\code{\link{MVoptrcd.maeT}}, \code{\link{MVoptrcd.maeA}}, \code{\link{Doptrcd.maeT}}, \code{\link{Doptrcd.maeA}}, 
\code{\link{Eoptrcd.maeT}} and \code{\link{Eoptrcd.maeA}}, respectively. Each function requires an initial connected row-column designs, 
generated using the function \code{\link{intcrcd.mae}}.  

The minimum value of \code{trt.N} and \code{col.N} is 3 and \code{trt.N} should be less than or equal to \code{col.N}. 
The linear fixed effects model results for given \code{trt.N} and \code{col.N} are obtained by setting \code{theta = 0.0}.

\code{Alg} specifies the array exchange and treatment exchange algorithm to be used that is adopted from Debusho, Gemechu and Haines (2016)   after adjusting for the row-column designs setup. If \code{Alg =} \code{"trtE"}, the function 
\code{optrcdmaeAT} perform  the treatment exchange procedure through deletion and addition of treatments at a time and selects a 
design with best treatment exchange with respect to the optimality criterion value. If \code{Alg =} \code{"arrayE"}, the function 
\code{optrcdmaeAT} perform the array exchange procedure through deletion and addition of candidate arrays at a time and selects a 
design with best array exchange with respect to the optimality criterion value.

\code{nrep} takes a value of greater than or equal to 2. However, to ensure optimality of the resultant design, 
the \code{nrep} should be greater than or equal to 10 and in addition, as \code{trt.N} and \code{col.N} increase, 
to ensure optimality of resultant design, it is advised to further increase the value of \code{nrep}
up to greater than or equal to 100. However, it has to be noted that as \code{trt.N} or \code{col.N} or
 \code{nrep} or all of them increase, computer time required to generate optimal or near-optimal
row-column design increases.

\code{itr.cvrgval} number of iterations during exchange procedure. It takes a value between 2 and \code{col.N}. It is used 
to speedup the computer search time by setting how long the user should wait for the exchange process to obtain any 
different (if any) design than the one that was produced as the result of the preceding exchange of the current array in the initial 
design with candidate array. This is mainly effective if \code{col.N} is very large. For example \code{itr.cvrgval = 2}, means the 
exchange procedure will jump to the next array test if the exchange of the two preceding arrays with candidate arrays results with the 
same efficient designs. The function  will not give error message if the users set \code{itr.cvrgval > col.N} and it will automatically 
set \code{itr.cvrgval = col.N}. The smaller the \code{itr.cvrgval} means the faster the exchange procedure is, but this will reduce the 
chance of getting optimal row-column design and users are advised to set \code{itr.cvrgval} closer to \code{col.N}. 


Remark: After the treatment exchange or array exchange procedure is completed, a dye-flip procedure is added to the internal functions of \code{optrcdmaeAT} stated above to further insure the optimality of the resulting optimal or near-optimal row-column designs. Thus, the procedure will flip (interchange) the treatments position within each array (column) and select the optimal dye-flip based on the optimality criteria of interest. This step is effective only for the large number of arrays and is efficient if \code{itr.cvgval} < \code{col.N} and there is a jump in the array exchange or treatment exchange procedure as stated above under the detail description of \code{itr.cvrgval}.  
}
\value{
Returns the resultant A-, MV-, D- or E-optimal or near-optimal row-column design with its corresponding score value and parametric combination 
saved in excel file in a working directory. In addition, the function \code{optrcdmaeAT} displays the graphical layout of the resultant 
optimal or near-optimal row-column designs. Specifically: 

\item{call}{the method call.}         
\item{v}{number of treatments.}
\item{b}{number of arrays.}
\item{theta}{theta value.}
\item{nrep}{number of replications of the optimization procedure.}  
\item{itr.cvrgval}{number of iterations required for convergence during the exchange procedure.}                          
\item{Optcrit}{optimality criteria.}                  
\item{Alg}{algorithm used.}
\item{OptdesF}{a \code{2 x col.N} obtained optimal or near-optimal row-column design.}
\item{Optcrtsv}{score value of the optimality criteria \code{'Optcrit'} of the resultant optimal or near-optimal row-column design \code{'OptdesF'}.}
\item{file_loc, file_loc2}{location where the summary of the resultant optimal or near-optimal row-column design is saved in .csv format.}
\item{equireplicate}{logical value indicating whether the resultant optimal or near-optimal row-column design is equireplicate or not.}
\item{vtrtrep}{vector of treatment replication of the resultant optimal or near-optimal row-column design.}     
\item{Cmat}{the C-matrix or  treatment information matrix of the  optimal or near-optimal row-column design.} 

The graphical layout of the resultant optimal or near-optimal row-column design.   

NB: The function \code{optrcdmaeAT} also saves the summary of the resultant optimal or near-optimal row-column design in .csv format in the working directory. 
Furthermore, the function reports only one final optimal or near-optimal row-column design, however, there is a possibility 
of more than one optimal or near-optimal row-column designs for a given parametric combination. 
The function \code{\link{graphoptrcd.mae}} can be used to view and rearrange the graphical layout of the resultant 
optimal or near-optimal row-column design on \code{tcltk} window. Alternative to the function \code{optrcdmaeAT}, a
GUI tcltk window can be used to generate optimal or near-optimal row-column designs, see \code{\link{mmenurcd.mae}} and \code{\link{fixparrcd.mae}}.   

}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal row-column designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2015). A-and D-optional row-column designs for two-colour cDNA microarray experiments using linear mixed effects models. \emph{South African Statistical Journal}, 49, 153-168.
}
\author{
Legesse Kassa Debusho, Dibaba Bayisa Gemechu, and Linda Haines
}
\seealso{
\code{\link{mmenurcd.mae}}, \code{\link{fixparrcd.mae}}, \code{\link{intcrcd.mae}}
}
\examples{
  \donttest{
  ##To obtain the A-optimal or near-optimal row-column design 
  ##using treatment exchange algorithm, set
  trt.N <- 3 #Number of treatments
  col.N <- 3 #Number of arrays
  theta <- 0 #theta value
  nrep <- 5  #Number of replications
  itr.cvrgval <- 6 #Number of iterations required during the exchange procedure
  Optcrit <- "A"   #Optimality criteria
  Alg <- "trtE"    #Algorithm
  
  Aoptrcdes <- optrcdmaeAT(trt.N = 3, col.N = 3, theta = 0, nrep = 5, 
                            itr.cvrgval = 6, Optcrit = "A", Alg = "trtE")
  
  summary(Aoptrcdes)
}
}

\keyword{A-optimal row-column designs}
\keyword{D-optimal row-column designs}
\keyword{E-optimal row-column designs}
\keyword{MV-optimal row-column designs}
\keyword{Microarray experiment} 
\keyword{Treatment exchange algorithm}
\keyword{Array exchange algorithm} 
