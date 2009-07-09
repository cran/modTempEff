\name{tempeff.fit}
\alias{tempeff.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of constrained segmented distributed lag model }
\description{
This is an internal function of package \code{modTempEff} and it should be not 
called by the user.
}
\usage{
tempeff.fit(y, X, Af = NULL, Ac = NULL, Xf = NULL, Xc = NULL, 
  V = NULL, ndx.seas, penalty = p.control(), only.seas = FALSE, 
  gam.fit.it = NULL, etastart = NULL, spstart = NULL, 
  fit.method = "magic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ See \code{\link{tempeff}} }
  \item{X}{ See \code{\link{tempeff}} }
  \item{Af}{ See \code{\link{tempeff}} }
  \item{Ac}{ See \code{\link{tempeff}} }
  \item{Xf}{ See \code{\link{tempeff}} }
  \item{Xc}{ See \code{\link{tempeff}} }
  \item{V}{ See \code{\link{tempeff}}}
  \item{ndx.seas}{ See \code{\link{tempeff}} }
  \item{penalty}{ See \code{\link{tempeff}} }
  \item{only.seas}{ See \code{\link{tempeff}} }
  \item{gam.fit.it}{ See \code{\link{tempeff}} }
  \item{etastart}{ See \code{\link{tempeff}} }
  \item{spstart}{ See \code{\link{tempeff}} }
  \item{fit.method}{ See \code{\link{tempeff}} }
}
\details{
This function is called by \code{tempeff} to fit the constrained segmented distributed lag model.
It is based on the function \code{gam.fit} of the \code{mgcv} package by S. Wood.
}
\value{
 A list of fit information.
}
\author{ Vito Muggeo }
\seealso{
\code{\link{tempeff}}
}

\keyword{regression}
