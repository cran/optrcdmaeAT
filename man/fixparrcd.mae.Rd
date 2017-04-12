\name{fixparrcd.mae}
\alias{fixparrcd.mae}
\title{
Sets parametric values
}
\description{
Creates a GUI \code{tcltk} window that allow the users to set or fix values for the parametric combinations to compute optimal or near-optimal row-column designs.}
\usage{
fixparrcd.mae(Optcrit)
}
\arguments{
  \item{Optcrit}{
character, specifying the optimality criteria to be used. \code{Optcrit} takes the letter \code{"A"}, \code{"MV"}, \code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal row-column designs, respectively.
}
}
\details{
\code{fixparrcd.mae} creates a pop-up GUI tcltk window that allow the users to set the parametric combinations to compute optimal or near-optimal row-column designs. The parameters include the number of treatments \code{trt.N}, number arrays  \code{col.N}, 
theta value  \code{theta}, number of replications of the optimization procedure \code{nrep} and number of iterations required during exchange procedure \code{itr.cvrgval}. 
Furthermore, on this window, the checkbox options that allow the users to choose whether or not they need to have the graphical 
layout of the resultant optimal or near-optimal row-column design, to make a choice between the two-alternative algorithms 
(treatment exchange and array exchange algorithms) and to print the summary of the resultant optimal or near-optimal row-column design on 
R console directly are available. 

After setting all the required parametric combinations and selecting the algorithm of interest, 
clicking on the search button on the set parametric combinations tcltk window, similar to the results that 
can be obtained when using the function \code{\link{optrcdmaeAT}}, the summary of the resultant optimal or near-optimal 
row-column design will be saved in the current working directory in .csv format and it will also be displayed on R console with 
graphical layout of the resultant optimal or near-optimal row-column designs.
}
\value{
The \code{fixparrcd.mae} function creates a pop-up tcltk window that allow the users to set 
the parametric combinations to compute optimal or near-optimal row-column designs.
}
\author{
Legesse Kassa Debusho, Dibaba Bayisa Gemechu, and Linda Haines
}
 
\seealso{
\code{\link{optrcdmaeAT}}, \code{\link{mmenurcd.mae}}, \code{\link{tcltk}}, \code{\link{TkWidgets}}
}
 