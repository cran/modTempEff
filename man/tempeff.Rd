\name{tempeff}
\alias{tempeff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Modelling temperature effects on mortality }
\description{
  Fits the constrained segmented distributed lag log-linear regression model to daily 
  time series data of mortality and temperature and additional confounding factors.
  }
\usage{
tempeff(formula, data, fcontrol = fit.control(), etastart = NULL, 
  drop.L, ...)
  }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ the model formula such as `response ~ parametric terms + csdl(temperature) + seas(day)', 
  see details.}
  \item{data}{ the dataset where the variables are stored. }
  \item{fcontrol}{ a list with components returned by \code{fit.control()}. }
  \item{etastart}{ possible starting values on the scale of the linear predictor. }
  \item{drop.L}{ integer, specifying whether the first `\code{drop.L}' observations should be removed
    before fitting. This is useful for model comparison purposes, see notes.}
  \item{\dots}{ additional arguments to be passed to \code{csdl()} in the formula, see details. }
}
\details{
  The function fits a log-linear regression model to assess the effects of temperature on mortality using a
  `constrained segmented distributed lag parameterization' (csdl). It is assumed that the data are daily time series 
  of mortality (or perhaps morbidity) and temperature. 
  The left hand side of the formula includes the response (daily counts), and the right hand side may include one
  or more of the following 
  \itemize{
  \item linear confounders (such as influenza epidemics or day-of-week);
  \item nonparametric long term trend, via the function \code{\link{seas}};
  \item the csdl effect of temperature via the function \code{\link{csdl}}.
  }
  All the arguments of \code{csdl()} may be passed via \code{...} directly in the call of \code{tempeff}. This may 
  facilitate the user when different models have to be fitted by changing only some of (and not all) the arguments
  of \code{csdl()}. See the example below.
  }
\value{
  The function returns an object of class \code{"modTempEff"}. It is the list returned by \code{gam.fit} of 
  package \code{mgcv} with the additional components
  \item{psi }{The estimated breakpoint with corresponding standard error (bayesian and frequentist).}
  \item{betaCold}{The estimated DL coefficients for the cold effect.}
  \item{SE.c}{The frequentist standard errors of the cold DL estimates. }
  \item{SE.c.bayes}{The bayesian standard errors of the cold DL estimates.}
  \item{ToTcold}{Estimate and frequentist standard error of the total (net) effect of cold.}
  \item{ToTcold.bayes}{Estimate and bayesian standard error of the total (net) effect of cold.}
  \item{edf.cold}{The df associated at each spline coefficient of the DL curve of cold.}
  \item{rank.cold}{The apparent dimension of the B-spline basis of the DL for cold.}
  \item{betaHeat}{The estimated DL coefficients for the heat effect.}
  \item{SE.h}{The frequentist standard errors of the heat DL estimates. }
  \item{SE.h.bayes}{The bayesian standard errors of the heat DL estimates.}
  \item{ToTheat}{Estimate and frequentist standard error of the total (net) effect of heat.}
  \item{ToTheat.bayes}{Estimate and bayesian standard error of the total (net) effect of heat.}
  \item{edf.heat}{The df associated at each spline coefficient of the DL curve of heat.}
  \item{rank.heat}{The apparent dimension of the B-spline basis of the DL for heat.}
  \item{rank.seas}{When \code{ndx.seas>0}, the apparent dimension of the B-spline basis for seasonality.}
  \item{edf.seas}{When \code{ndx.seas>0}, the df associated at spline coefficients of seasonality.}
  \item{fit.seas}{When \code{ndx.seas>0}, the fitted long-term trend (on the log scale).}
}
\references{
Muggeo, V.M.R. (2008) Modeling temperature effects on mortality: multiple 
    segmented relationships with common break points
    \emph{Biostatistics} \bold{9}, 613--620.

Muggeo, V.M.R. (2009) Analyzing Temperature Effects on Mortality Within the R 
  Environment: The Constrained Segmented Distributed Lag Parameterization
    \emph{Journal of Statistical Software}, \bold{32} 12, 1--17.
    }

\author{Vito Muggeo, \email{vito.muggeo@unipa.it}}

\note{
  When a \code{csdl} term is included in the formula, the first \code{max(L)} observations are discarded 
  before model fitting. When a csdl term is \emph{not} included, the argument \code{drop.L} may be used to discard 
  the first \code{drop.L} observations anyway. Fitting models with the same number of observations may be useful 
  to compare them via likelihood-based criteria (via \code{anova.modTempEff}, say). \cr
  \code{tempeff()} returns objects of class \code{"modTempEff"}, so proper methods may be employed. The returned 
  object has class \code{"modTempEff"} even if \code{tempeff()} is called without \code{csdl()} in the formula, 
  or even if the model is fitted with fixed (not estimated) breakpoints (via \cr 
  \code{tempeff(..,fcontrol=fit.control(it.max=0))}). 
  }

\seealso{\code{\link{modTempEff-package}}, \code{\link{plot.modTempEff}}, \code{\link{summary.modTempEff}}, 
\code{\link[mgcv]{gam.fit}} in package \code{mgcv}}

%\seealso{\code{\link{ivreg.fit}}, \code{\link[stats]{lm}}, \code{\link[stats]{lm.fit}}}

\examples{
\dontrun{
library(modTempEff)
data(dataDeathTemp)
o1<-tempeff(dec1~day+factor(dweek)+factor(year)+factor(month)+
      csdl(mtemp,L=c(60,60),psi=20),
      data=dataDeathTemp, fcontrol = fit.control(display=TRUE))

#add a ridge penalty: note how you *can* specify ridge!
#you do NOT need to use csdl(..,ridge=..) 
o2<-update(o1, ridge=list(cold="l^2", heat="l^2"))

#a model without temperature effects (the first drop.L obs are dropped)
o3<-tempeff(dec1~day+factor(dweek)+factor(year)+factor(month),
      data=dataset,drop.L=60)

#see ?anova.modTempEff for model comparisons
}
}

\keyword{regression}