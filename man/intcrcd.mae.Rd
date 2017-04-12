\name{intcrcd.mae}
\alias{intcrcd.mae}
\title{
Generates initial connected row-column design
}
\description{
Generates a random initial connected row-column design for a given number of arrays \code{b} of size \code{k = 2} and the number of treatments \code{v}.}
\usage{
intcrcd.mae(trt.N, col.N)
}
\arguments{
 \item{trt.N}{
integer, specifying number of treatments, \code{v}. 
}
  \item{col.N}{
integer, specifying number of arrays, \code{b}.
}
}
\value{
Returns a \code{2 x b} connected row-column design with \code{b} arrays of size \code{k = 2} and number of treatments \code{v}. 
}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.

Gemechu, D. B., Debusho, L. K., and Haines, L. M. (2015). A-and D-optional row-column designs for two-colour cDNA microarray experiments using linear mixed effects models. \emph{South African Statistical Journal}, 49, 153-168.
}
\author{
Legesse Kassa Debusho, Dibaba Bayisa Gemechu, and Linda Haines
}

\seealso{
\code{\link{optrcdmaeAT}}, \code{\link{cmatrcd.mae}}
}
\examples{
    #Initial connected row-column design for 
    trt.N <- 4 #Number of treatments
    col.N <- 4 #Number of arrays

    intcrcd.mae(trt.N = 4, col.N = 4)
}
\keyword{Initial row-column design}
\keyword{Connected row-column design}% __ONLY ONE__ keyword per line
