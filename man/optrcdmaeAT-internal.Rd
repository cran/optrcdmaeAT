\name{optrcdmaeAT-internal}
\alias{optrcdmaeAT-internal}
\alias{Aoptrcd.maeA}
\alias{Aoptrcd.maeT}
\alias{MVoptrcd.maeA}
\alias{MVoptrcd.maeT}
\alias{Doptrcd.maeA}
\alias{Doptrcd.maeT}
\alias{Eoptrcd.maeA}
\alias{Eoptrcd.maeT}
\title{
Internal functions
}
\description{
Functions for internal usage only.
}
\usage{
 
## Computes A-optimal or near-optimal row-column designs
## using array exchange algorithm
Aoptrcd.maeA(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes A-optimal or near-optimal row-column designs
## using treatment exchange algorithm
Aoptrcd.maeT(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes MV-optimal or near-optimal row-column designs
## using array exchange algorithm
MVoptrcd.maeA(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes MV-optimal or near-optimal row-column designs
## using treatment exchange algorithm
MVoptrcd.maeT(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes D-optimal or near-optimal row-column designs
## using array exchange algorithm
Doptrcd.maeA(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes D-optimal or near-optimal row-column designs
## using treatment exchange algorithm
Doptrcd.maeT(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes E-optimal or near-optimal row-column designs
## using array exchange algorithm
Eoptrcd.maeA(trt.N, col.N, theta, nrep, itr.cvrgval)
 
## Computes E-optimal or near-optimal row-column designs
## using treatment exchange algorithm
Eoptrcd.maeT(trt.N, col.N, theta, nrep, itr.cvrgval)
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
integer, specifying number of iterations required for convergence during the exchange procedure. See \code{\link{optrcdmaeAT}} documentation for details.
}
}
\details{
These functions are handled via a generic function \code{\link{optrcdmaeAT}}. Please refer to the \code{\link{optrcdmaeAT}} documentation for details. 
}

\references{
Debusho, L. K., Gemechu, D. B. and Haines, L. (2018). Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Communications in Statistics - Simulation and Computation, https://doi.org/10.1080/03610918.2018.1429617}.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2015). A-and D-optional row-column designs for two-colour cDNA microarray experiments using linear mixed effects models. \emph{South African Statistical Journal}, 49, 153-168.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\seealso{
\code{\link{optrcdmaeAT}}
}
\keyword{A-optimal row-column designs}
\keyword{D-optimal row-column designs}
\keyword{E-optimal row-column designs}
\keyword{MV-optimal row-column designs}
\keyword{Microarray experiment} 
\keyword{Treatment exchange algorithm}
\keyword{Array exchange algorithm} 
