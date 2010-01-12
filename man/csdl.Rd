\name{csdl}
\alias{csdl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Defining the constrained segmented distributed lag term
}
\description{
Function employed within the \code{\link{tempeff}} formula to specify the variable with a csdl relationship. 
The function does not fit the model, it simply returns information exploited by \code{\link{tempeff}} to fit a
\code{"modTempEff"} model.
}
\usage{
csdl(z, psi, L, ridge = NULL, ndx = round(L/3),
    DL = FALSE, diff.varying = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{the variable, typically the temperature, having the constrained segmented distributed
  lag parameterization with the response.
  }
  \item{psi}{numeric to provide the starting value for one or two breakpoints of the 
    constrained segmented distributed lag relationship.
  }
  \item{L}{a numerical two-length vector to specify how many lags have to be considered to assess
  the effect of cold and heat.
  }
  \item{ndx}{\emph{apparent} dimension (i.e. the rank) of the two B-spline bases for the DL curves 
  of cold and heat. Default to \code{ndx=round(L/3)}.
  }
  \item{ridge}{a two-length named list of characters to specify the possible ridge penalty to be applied
  to DL coefficients. This list has to be named ('cold' and 'heat') and each component has to be
  in terms of \code{"l"}. Use \code{ridge=NULL} whether no ridge penalty has to be employed. See examples below.
  }
  \item{DL}{logical indicating if the difference penalty should be applied to the DL coefficients
    or to spline coefficients. Default to spline coefficients (\code{DL=FALSE}).
  }
  \item{diff.varying}{logical indicating if the difference penalty should be global 
  or depending on the lag value to penalize differences mainly at larger rather than early lags.
  }
}
\details{
This function has to be used within the formula in \code{\link{tempeff}}. It returns information to 
fit a constrained segmented distributed lag parameterization within the Poisson regression model 
via the function \code{\link{tempeff}}.
}
\value{
A list with the temperature variable and the arguments as attributes.
}
\note{
All the arguments of \code{csdl} may be passed via \code{...} in the call of \code{tempeff}. This feature may 
be useful when the same model has to be re-fitted (via \code{update}) by modifying only one argument 
of \code{csdl()}. See \code{\link{tempeff}}.
}
\references{ 
Eilers, P., Marx, B. D. (1996). Flexible Smoothing with B-splines and Penalties 
\emph{Statistical Science} \bold{11}, 89--121.

Wood, S. N. (2006) \emph{Generalized Additive Models: An Introduction with R.} Chapman and Hall/CRC Press.
                                     
Muggeo, V. M. R. (2008) Modeling temperature effects on mortality: multiple 
segmented relationships with common break points
\emph{Biostatistics} \bold{9}, 613--620.
 }
\author{
Vito Muggeo
}
\seealso{\code{\link{tempeff}}}
\examples{
\dontrun{
# Evaluate temperature effects up to 45 lags for cold and heat and 
# a single breakpoint; use 20 as starting value, a global difference 
# penalty on spline coefficients and no ridge penalty
csdl(my.temperature,psi=20,L=c(45,45),ridge=NULL)

# Evaluate temperature effects up to 45 lags for cold and 15 lag for 
# heat, via P-splines with a global difference penalty on DL coefficients 
# and an additional quadratic ridge penalty
csdl(my.temperature,psi=20,L=c(45,15),DL=TRUE,ridge=list(cold="l^2", 
  heat="l^2"))
}
}
\keyword{ regression }
