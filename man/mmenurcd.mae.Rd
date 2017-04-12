\name{mmenurcd.mae}
\alias{mmenurcd.mae}
\title{
Creates the main menu tcltk window
}
\description{
Alternative to directly using the function \code{\link{optrcdmaeAT}} on the R console, the function \code{mmenurcd.mae} creates the main menu tcltk window with the list of four optimality criteria (A, MV, D and E) for search of optimal or near-optimal row-column designs for two-colour cDNA microarray experiments. }
\usage{
mmenurcd.mae()
}
\details{
The main menu window created using the function \code{mmenurcd.mae} contains the list of four optimality criteria (the A-,
MV -, D- and E) for row-column designs that allows user to choose an optimality criterion of interest. The function 
will then call for a function \code{fixparrcd.mae} for further option regarding the setup of parametric combinations and output, see \code{fixparrcd.mae} documentation. 
}
\author{
Legesse Kassa Debusho, Dibaba Bayisa Gemechu, and Linda Haines
}

\seealso{
\code{\link{optrcdmaeAT}}, \code{\link{fixparrcd.mae}}, \code{\link{tcltk}}, \code{\link{TkWidgets}}
}
