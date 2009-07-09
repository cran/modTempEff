\name{plot.modTempEff}
\alias{plot.modTempEff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for the class 'modTempEff'}
\description{
  Plots distributed lags curves from the \code{modTempEff} fit.
}
\usage{
\method{plot}{modTempEff}(x, which = c("cold", "heat"), new = TRUE, 
  var.bayes = FALSE, add = FALSE, delta.rr = TRUE, level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{"modTempEff"}. }
  \item{which}{ Which DL curve should be plotted? for cold, heat or both of them (default). }
  \item{new}{ logical indicating if a new device should be opened.}
  \item{var.bayes}{ logical indicating if the 'Bayesian' rather than the frequentist
    standard errors should be employed to compute the pointwise confidence intervals to be plotted}
  \item{add}{ an alias for \code{new}. }
  \item{delta.rr}{ logical indicating if the DL curves should be plotted on the log scale or on 
  the per cent change in relative risk, i.e. 100*(exp(.)-1).}
  \item{level}{ the selected confidence level of the pointwise confidence intervals to be plotted}
  \item{\dots}{ additional arguments.. }
}
\details{
     Takes a fitted 'modTempEff' object produced by 'tempeff()' and plots the 
     DL curves for cold and heat effect with relevant pointwise confidence intervals. 
}
\value{
The function simply plots the required estimated DL curve.
}
%\references{%}
\author{ Vito Muggeo }
\seealso{ \code{\link{tempeff}}}
\examples{
  \dontrun{
  #obj is an object returned by tempeff()
  #plots DL curves for cold and heat with 95% pointwise confidence 
  # using frequentist standard errors
  plot(obj)
  
  #plots the estimated DL curve only for heat with 90% pointwise confidence 
  # using bayesian standard errors
  plot(obj, "heat", var.bayes=TRUE, level=.90)
  }
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
