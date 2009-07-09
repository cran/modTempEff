\name{tempeff}
\alias{tempeff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Modelling Temperature Effects on Mortality }
\description{
  Fits the constrained segmented distributed lags log-linear regression model to daily 
  time series data of mortality and temperature and additional confounding factors.
  }
\usage{
tempeff(formula, z, data, tcontrol = temp.control(), 
  pcontrol = p.control(), fcontrol = fit.control(), 
  etastart = NULL, ndx.seas = 0, ...)
  }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ the model formula including the `response ~ parametric_terms', see details.}
  \item{z}{ the temperature variable. }
  \item{data}{ the dataset where the variables are stored. }
  \item{tcontrol}{ a list with components returned by \code{temp.control()}. }
  \item{pcontrol}{ a list with components returned by \code{p.control()}. }
  \item{fcontrol}{ a list with components returned by \code{fit.control()}. }
  \item{etastart}{ possible starting values on the scale of the linear predictor. }
  \item{ndx.seas}{ possible apparent dimension of the B-spline basis for seasonality 
  (actually the effective dimension is ndx.seas+3). If \code{ndx.seas=0} no spline for seasonality is assumed. }
  \item{\dots}{ additional arguments to be passed to \code{tempeff.fit()}; currently unimplemented. }
}
\details{
  This function fits a log-linear regression model to assess the effects of temperature on mortality.
  It is assumed that the data are daily time series of mortality (or perhaps morbidity) and temperature. 
  The response and the confounders (such as influenza epidemics or day-of-week) have to be specified in the 
  \code{formula} and the temperature itself in specified in the argument \code{z}. Long-term trend and 
  seasonality may be modelled via P-splines by specifying the apparent dimension of the relevant basis via 
  the argument \code{ndx.seas}.
  }
\value{
  The function returns an object of class \code{"modTempEff"}. It is the list returned by \code{gam.fit} of 
  package \code{mgcv} with the additional components
  \item{psi }{The estimated breakpoint with corresponding standard error (bayesian and frequentist).}
  \item{betaFreddo}{The estimated DL coefficients for the cold effect.}
  \item{SE.f}{The frequentist standard errors of the cold DL estimates. }
  \item{SE.f.bayes}{The bayesian standard errors of the cold DL estimates.}
  \item{ToTfreddo}{Estimate and frequentist standard error of the total (net) effect of cold.}
  \item{ToTfreddo.bayes}{Estimate and bayesian standard error of the total (net) effect of cold.}
  \item{edf.cold}{The df associated at each spline coefficient of the DL curve of cold.}
  \item{rank.cold}{The apparent dimension of the B-spline basis of the DL for cold.}
  \item{betaCaldo}{The estimated DL coefficients for the heat effect.}
  \item{SE.c}{The frequentist standard errors of the heat DL estimates. }
  \item{SE.c.bayes}{The bayesian standard errors of the heat DL estimates.}
  \item{ToTcaldo}{Estimate and frequentist standard error of the total (net) effect of heat.}
  \item{ToTcaldo.bayes}{Estimate and bayesian standard error of the total (net) effect of heat.}
  \item{edf.heat}{The df associated at each spline coefficient of the DL curve of heat.}
  \item{rank.heat}{The apparent dimension of the B-spline basis of the DL for heat.}
  \item{rank.seas}{When \code{ndx.seas>0}, the apparent dimension of the B-spline basis for seasonality.}
  \item{edf.seas}{When \code{ndx.seas>0}, the df associated at spline coefficients of seasonality.}
  \item{fit.seas}{When \code{ndx.seas>0}, the fitted long-term trend (on the log scale).}
}
\references{
Muggeo, V.M.R. (2008) Modeling temperature effects on mortality: multiple 
    segmented relationships with common break points
    \emph{Biostatistics} \bold{9}, 613--620.}
\author{ Vito Muggeo, \email{vito.muggeo@unipa.it}}
\note{ The first 'max(L)' observations are discarded before model fitting.}

\seealso{\code{\link{modTempEff-package}}, \code{\link{plot.modTempEff}}, \code{\link{summary.modTempEff}}, 
\code{\link{gam.fit}} in package \code{mgcv}
}

\examples{
\dontrun{
library(modTempEff)
data(dataset)
o1<-tempeff(dec1~day+factor(dweek)+factor(year)+factor(month), data=dataset,
      tcontrol = temp.control(psi=20, L=c(60,60)), z=mtemp,
      pcontrol = p.control(ridge.formulas=NA),
      fcontrol = fit.control(display=TRUE))

o2<-update(o1A,
  pcontrol=p.control(ridge.formulas=list(cold="xlag^3", heat="xlag^3")))
}
}