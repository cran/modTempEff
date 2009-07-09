\name{temp.control}
\alias{temp.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary function to control the fit of the temperature effect}
\description{
Specifies some inputs to fit the segmented effect of temperature, including starting values
for the threshold value, the maximum lags for cold and heat, and the `apparent' dimension of
the B-spline bases for cold and heat.
}
\usage{
temp.control(psi, L = c(60, 60), ndx = round(L/3), heat.power = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{Numeric to provide the starting value of the one or two breakpoint of the segmented relationship.}
  \item{L}{ A numerical two-length vector to specify how many lags have to be considered to assess
  the effect of cold and heat.}
  \item{ndx}{ The apparent dimension of the two B-spline bases for the DL curves of cold and heat. }
  \item{heat.power}{ Currently not implemented. }
}
\details{
This auxiliary function is employed to specify the starting value for the
one or two breakpoints, the maximum lags for cold and heat effect, and the `apparent' dimension,
i.e. the number of columns, of the B-spline bases to fit the DL curves of cold and heat.
}
\value{
  A list with the arguments as components to be used by 'tempeff' or 'tempeff.fit'.
}
\references{ 
Muggeo, V.M.R. (2008) Modeling temperature effects on mortality: multiple 
    segmented relationships with common break points
    \emph{Biostatistics} \bold{9}, 613--620.
}
\author{ Vito Muggeo }
\seealso{ \code{\link{tempeff}}}
\examples{
\dontrun{
# Evaluate temperature effects up to 45 lags and a single breakpoint; 
# use 20 as starting value
temp.control(psi=20,L=c(45,45))
}
}

\keyword{regression}
