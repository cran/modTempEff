\name{modTempEff-package}
\alias{modTempEff-package}
\alias{modTempEff}
\docType{package}
\title{
  Modelling temperature effects on mortality via constrained segmented 
  distributed lag parameterization}
\description{
  Estimation of a Poisson log linear regression to model the temperature effect on health 
  using the `constrained segmented distributed lag parameterization' which allows to account simultaneously 
  for typical features of temperature effects: nonlinear effect delayed over several days.
}
\details{
\tabular{ll}{
Package: \tab modTempEff\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-07-08\cr
License: \tab GPL\cr
}
Package \code{modTempEff} fits a Poisson log linear regression to model the temperature effect on health 
  using the `constrained segmented distributed lag parameterization' which allows to account simultaneously 
  for typical features of temperature effects: nonlinear effect delayed over several days. 
  Fitting is performed by the function \code{tempeff} and print, summary, and 
  plot methods are included to summarize results.
}
\author{
Vito M.R. Muggeo \email{vito.muggeo@unipa.it}
%##Maintainer: Vito M.R. Muggeo <vito.muggeo@unipa.it>
}
\references{
Muggeo, V.M.R. (2008) Modeling temperature effects on mortality: multiple 
    segmented relationships with common break points
    \emph{Biostatistics} \bold{9}, 613--620.

Muggeo, V.M.R., Hajat, S. (2009) Modelling the nonlinear multiple-lag effects of 
  ambient temperature on mortality in Santiago and Palermo : a constrained 
  segmented distributed lag approach  \emph{Occupational Environmental Medicine} 
  to appear.

Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points.
    \emph{Statistics in Medicine} \bold{22}, 3055--3071.


Some background references on temperature effect on mortality via time series data (epidemiological papers)

Armstrong, B. (2006) Models for the relationship between ambient temperature and daily
mortality \emph{Epidemiology} \bold{17}, 624--631.

Basu, R., Samet, J. (2002). Relation between elevated ambient temperature and mortality: a
review of the epidemiologic evidence  \emph{Epidemiological Reviews} \bold{24}, 190--202.
}

\keyword{ package }
\keyword{ models }
\seealso{
\code{\link[<mgcv>:<mgcv>-package]{mgcv}}
}
\examples{
\dontrun{
data(dataset)
o1<-tempeff(dec1~day+factor(dweek)+factor(year)+factor(month), data=dataset,
      tcontrol = temp.control(psi=20, L=c(60,60)), z=mtemp,
      pcontrol = p.control(ridge.formulas=list(cold="xlag^2",heat="xlag^2")),
      fcontrol = fit.control(display=TRUE))
}
}
