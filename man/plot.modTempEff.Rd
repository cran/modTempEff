\name{plot.modTempEff}
\alias{plot.modTempEff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for the class 'modTempEff'}
\description{
  Plots distributed lags curves from the \code{modTempEff} fit.
}
\usage{
\method{plot}{modTempEff}(x, which = c("cold", "heat"), add=FALSE, new=TRUE, 
  var.bayes = FALSE, delta.rr = TRUE, level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{"modTempEff"}. }
  \item{which}{ Which DL curve should be plotted? for cold, heat or both of them (default). }
  \item{add}{ logical; if \code{TRUE} the fitted DL curve for cold \emph{or} heat is added 
    to an existing plot. }
  \item{new}{ logical indicating if a new device should be opened. If \code{add=TRUE}, \code{new}
    is set to \code{FALSE}.}
  \item{var.bayes}{ logical indicating if the 'Bayesian' rather than the frequentist
    standard errors should be employed to compute the pointwise confidence intervals to be plotted}
  \item{delta.rr}{ logical indicating if the DL curves should be plotted on the log scale or as
    per cent change in relative risk, i.e. 100*(exp(.)-1).}
  \item{level}{ the selected confidence level of the pointwise confidence intervals to be plotted}
  \item{\dots}{ additional arguments.. }
}
\details{
     Takes a fitted \code{"modTempEff"} object produced by \code{tempeff()} and plots the 
     DL curves for cold and heat effect with relevant pointwise confidence intervals.
     \code{plot.modTempEff} also works with objects with fixed (not estimated) breakpoint, namely 
     fits returned by \cr 
     
     \code{tempeff(.., fcontrol=fit.control(it.max=0))}. \cr 
     
     Note \code{add=TRUE} makes sense (and works) only for a single (cold \emph{or} heat) DL curve to be superimposed to an existing plot. 
}
\value{
The function simply plots the required estimated DL curve. If the fitted model includes only a smooth term for 
the long term trend, \code{plot.modTempEff} draws it.
}
%\references{%}
\author{ Vito Muggeo }
\seealso{ \code{\link{tempeff}}}
\examples{
  \dontrun{
  #obj is an object returned by tempeff()
  #plots DL curves for cold and heat with 95\% pointwise CI
  # using frequentist standard errors
  plot(obj)
  
  #plots the estimated DL curve only for heat with 90\% pointwise CI
  # using bayesian standard errors
  plot(obj, "heat", var.bayes=TRUE, level=.90)
  }
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
