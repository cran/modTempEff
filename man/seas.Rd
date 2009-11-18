\name{seas}
\alias{seas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specifying a smooth long term trend within a `modTempEff' model 
}
\description{
Function employed within the \code{\link{tempeff}} formula to specify the variable with a smooth effect. 
The function does not fit the model, it simply returns information exploited by \code{\link{tempeff}}.
}
\usage{
seas(x, ndx = stop("please, provide `ndx' in seas()"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the long term trend variable
  }
  \item{ndx}{the apparent dimension (i.e. the rank) of the B-spline employed to 
  model the long trend.
  }
}
\details{
The function is used within the formula of \code{tempeff()} when the long term trend of the daily time series 
is fitted via a smooth term. A standard P-spline with a cubic B-spline basis and a second-order 
difference penalty is employed.
}
\value{
A list with relevant information.
}
\references{
Eilers, P., Marx, B. (1996). Flexible Smoothing with B-splines and Penalties 
\emph{Statistical Science} \bold{11}, 89--121.

Wood, S.N. (2006) \emph{Generalized Additive Models: An Introduction with R.} 
Chapman and Hall/CRC Press.
}
\author{
Vito Muggeo
}
\seealso{\code{\link{csdl}}, \code{\link{tempeff}}}
%\examples{
%}
\keyword{ regression }
