\name{p.control}
\alias{p.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary function to control the P-spline fit for distributed lags}
\description{
  Auxiliary function to detail the options concerning the P-splines employed to 
  smooth the distributed lags (DL) curves of cold and heat.
}
\usage{
p.control(DL = FALSE, diff.varying = FALSE, 
  ridge.formulas = list(cold = "xlag^2", heat = "xlag^2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DL}{logical indicating if the difference penalty should be applied to the DL coefficients
    or spline coefficients. Default to spline coefficients (\code{DL=FALSE})}
  \item{diff.varying}{ logical indicating if the difference penalty should be global 
  or depending on the lag, i.e. to penalize differences mainly at longer lags. }
  \item{ridge.formulas}{ A two-length named list of characters to specify the possible ridge penalty to be applied
  to DL coefficients. This list has to be named ('cold' and 'heat') and each component has to be
  in terms of \code{"xlag"}. See examples below.}
}
\details{
  DL curves are estimated via P-splines, namely B-spline basis with some penalty on coefficients 
  to guarantee smoothness. The difference penalty ensures smooth curve over the lag range and 
  the possible ridge penalty allows the DL curve approaching to zero at longer lags.   
}
\value{
  A list with the arguments as components to be used by 'tempeff' or 'tempeff.fit'.
}
\references{ 
Eilers, P., Marx, B. (1996). Flexible Smoothing with B-splines and Penalties 
\emph{Statistical Science} \bold{11}, 89--121.
                                     
Muggeo, V.M.R. (2008) Modeling temperature effects on mortality: multiple 
segmented relationships with common break points
\emph{Biostatistics} \bold{9}, 613--620.
 }
\author{Vito Muggeo}
\seealso{ 
\code{\link{tempeff}}
}
\examples{
\dontrun{
#P-spline fit with a global difference penalty on DL coefficients and 
#an additional quadratic ridge penalty
p.control(DL=FALSE,diff.varying=FALSE,
  ridge.formulas=list(cold="xlag^2",heat="xlag^2"))

#P-spline fit with a global difference penalty of the DL coefficients 
#without ridge penalty
p.control(DL=TRUE,diff.varying=FALSE,ridge.formulas=NA)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
