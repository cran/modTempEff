\name{coef.modTempEff}
\alias{coef.modTempEff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract DL coefficients from a `modTempEff' object
}
\description{
The function extracts the estimated coefficients of the DL curves for cold and/or
heat effect.
}
\usage{
coef.modTempEff(object, which = c("cold", "heat"), L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
the \code{"modTempEff"} object returned by \code{\link{tempeff}}.
}
  \item{which}{
which DL curve should be returned? 
}
  \item{L}{
the number of DL coefficients required. `\code{L}+1' coefficients for lags 0 to \code{L} are returned.
}
  \item{\dots}{
additional arguments (ignored).
}
}
\details{
The resulting estimates, returned as a matrix, are DL coefficients for the cold and/or the heat effect. Each coefficient
  at specific lag represents the log relative risk (of mortality) for one-unit increase in cold (or heat) values.
  }
\value{
A matrix with DL coefficients. Each row corresponds to a specific lag.
}
\author{
Vito Muggeo
}
\seealso{
\code{\link{tempeff}}
}
\examples{
\dontrun{
 #continues from ?tempeff
 coef(o1,"heat",L=7) #log RR for heat corresponding to lag 0 to 7
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
