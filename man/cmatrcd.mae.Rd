\name{cmatrcd.mae}
\alias{cmatrcd.mae}
\title{
Computes the treatment information matrix 
}
\description{
Computes the information matrix (C-matrix) for treatment effects under either the linear fixed effects model  or the linear mixed effects model setting for a given row-column design.
}
\usage{
cmatrcd.mae(trt.N, col.N, theta, des)
}
\arguments{
  \item{trt.N}{
integer, specifying number of treatments, \code{v}. 
}
  \item{col.N}{
integer, specifying number of arrays (columns), \code{b}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{des}{
matrix, a \code{2 x b} row-column design with \code{b} arrays/columns of size \code{k = 2} and \code{v} treatments.
}
}
\value{
Returns a \code{v x v} treatment information matrix (C-matrix).}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2015). A-and D-optional row-column designs for two-colour cDNA microarray experiments using linear mixed effects models. \emph{South African Statistical Journal}, 49, 153-168.
}
\author{
Legesse Kassa Debusho, Dibaba Bayisa Gemechu, and Linda Haines
}

\seealso{
\code{\link{optrcdmaeAT}}, \code{\link{fixparrcd.mae}}, \code{\link{intcrcd.mae}}
}
\examples{

##Information matrix

     trt.N <- 3 
     col.N <- 3 
     theta <- 0.1 
     rcdes <- intcrcd.mae(trt.N = 3, col.N = 3)

     cmatrcd.mae(trt.N = 3, col.N = 3, theta = 0.1, des = rcdes)
}
\keyword{Information matrix}
\keyword{C-matrix} 
