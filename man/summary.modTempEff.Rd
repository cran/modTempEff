\name{summary.modTempEff}
\alias{summary.modTempEff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for the class 'modTempEff'}
\description{
Summarizes fit for the constrained segmented distributed effect model.}
\usage{
\method{summary}{modTempEff}(object, spar = TRUE, 
  digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class \code{"modTempEff"} }
  \item{spar}{ logical indicating if values of smoothing parameters should be printed }
  \item{digits}{ number of digits to be printed }
  \item{\dots}{ additional arguments.. }
}
\details{
  Prints the most important features of a \code{modTempEff} object including fit summary 
  (AIC, BIC, Ubre, residual deviance) and point estimates along with standard errors of the net effect
  of cold and heat, and the breakpoints where mortality reaches its minimum. Smoothing parameters 
  (selected via \code{gam.fit()} of the mgcv package) are also printed when \code{spar=TRUE}.
}
\author{ Vito Muggeo}
\seealso{ \code{\link{print.modTempEff}} }
\keyword{ regression}
