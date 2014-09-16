\name{logLik.modTempEff}
\alias{logLik.modTempEff}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract the log likelihood for a modTempEff fit}
\description{Function to extract the log-likelihood for a fitted \code{modTempEff}
model.
}
\usage{
\method{logLik}{modTempEff}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{object}{ fitted model objects of class \code{modTempEff} as produced by \code{tempeff()}.}
\item{...}{ignored}

}
\details{ This function extracts the \emph{penalized} log likelihood for a fit produced by \code{tempeff()},
namely a Poisson GAM fit as returned by \code{mgcv::gam()}.

Notice that the model degrees of freedom are the effective degrees of freedom
and not the number of coefficients, as the model is estimated by penalized likelihood.

}

\value{ Standard \code{logLik} object: see \code{\link{logLik}}.
}


%\references{Hastie and Tibshirani, 1990, Generalized Additive Models.
%Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for
%generalized additive models. J.R.Statist. Soc. B 70(3):495-518
%}

\author{ Vito M. R. Muggeo \email{vito.muggeo@unipa.it} based directly on \code{logLik.gam} in \code{mgcv} 
  package by S. Wood.}


\seealso{  \code{\link{AIC}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






