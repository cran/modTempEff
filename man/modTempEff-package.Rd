\name{modTempEff-package}
\alias{modTempEff-package}
\alias{modTempEff}
\docType{package}
\title{
  Modelling temperature effects on mortality via the constrained segmented 
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
Version: \tab 1.5\cr
Date: \tab 2010-01-12\cr
License: \tab GPL\cr
}
Package \code{modTempEff} fits a Poisson log linear regression to model the temperature effect on health 
  using the `constrained segmented distributed lag parameterization' which allows to account simultaneously 
  for typical features of temperature effects: nonlinear effect delayed over several days. 
  Fitting is performed by the function \code{tempeff} and 'print', 'summary', 'plot', 'anova' 
  and 'coef' methods are included to summarize and to display results.
}
\author{
Vito M.R. Muggeo \email{vito.muggeo@unipa.it}
%##Maintainer: Vito M.R. Muggeo <vito.muggeo@unipa.it>
}
\references{
Muggeo, V.M.R. (2010) Analyzing temperature effects on mortality within the R 
  environment: the constrained segmented distributed lag parameterization
    \emph{Journal of Statistical Software} \bold{32}, 12, 1--17.

Muggeo, V.M.R. (2008) Modeling temperature effects on mortality: multiple 
    segmented relationships with common break points
    \emph{Biostatistics} \bold{9}, 613--620.

Muggeo, V.M.R., Hajat, S. (2009) Modelling the nonlinear multiple-lag effects of 
  ambient temperature on mortality in Santiago and Palermo : a constrained 
  segmented distributed lag approach. \emph{Occupational Environmental Medicine} 
  \bold{66}, 584--591. 

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
\code{\link[mgcv]{mgcv}}
}
\examples{
\dontrun{
data(dataDeathTemp)
o1<-tempeff(dec1~day+factor(dweek)+factor(year)+factor(month)+
      csdl(mtemp,L=c(60,60),psi=20, ridge=list(cold="l^2",heat="l^2")),
      data=dataDeathTemp, fcontrol = fit.control(display=TRUE))

o2<-tempeff(dec1~seas(day,30)+
      csdl(mtemp,L=c(60,60),psi=20, ridge=list(cold="l^2",heat="l^2")),
      data=dataDeathTemp, fcontrol = fit.control(display=FALSE))

}
}
