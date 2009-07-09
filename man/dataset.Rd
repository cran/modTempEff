\name{dataset}
\alias{dataset}
\docType{data}
\title{Simulated dataset }
\description{
  Simulated dataset of daily time series of mortality, temperature, pollutant, 
  and seasonal variables
}
\usage{data(dataset)}
\format{
  A data frame with 1825 observations on the following 8 variables.
  \describe{
    \item{\code{dec1}}{daily counts of mortality}
    \item{\code{mtemp}}{daily mean temperatute}
    \item{\code{month}}{month}
    \item{\code{year}}{year}
    \item{\code{day}}{day}
    \item{\code{dweek}}{day of week}
    \item{\code{decNS}}{another daily counts of mortality}
    \item{\code{dec2}}{yet another daily counts of mortality}
  }
}
\details{
  These data represent a typical dataset employed to investigate the short term 
  effect of temperature (and pollutant) on health via daily time series.
}
%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
\dontrun{data(dataset)}
}
\keyword{datasets}
