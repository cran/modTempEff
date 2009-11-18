\name{anova.modTempEff}
\alias{anova.modTempEff}
\title{The anova method for a 'modTempEff' object}
\description{
Comparing \code{"modTempEff"} objects returned by \code{tempeff()} using an analysis of deviance table.
}
\usage{
\method{anova}{modTempEff}(object, ..., dispersion = NULL, test = NULL)
}
\arguments{
  \item{object,...}{ fitted model objects of class \code{"modTempEff"} returned by \code{tempeff()}.}
  \item{dispersion}{ currently ignored. }
  \item{test}{ what sort of likelihood based criterion has to be used for model comparisons. One of
          \code{"Chisq"}, \code{"F"} or \code{"Cp"}.}
%  \item{\dots}{ additional arguments.. }
}
\details{\code{anova.modTempEff} performs model comparisons in terms of likelihood-based criteria 
depending on its argument \code{test}. When \code{test="Chisq"}, likelihood ratio tests are carried out; 
however note that the p-values for smooth terms (included via \code{seas()} or \code{csdl()}) are approximate only.
The function does \emph{not} work for a single \code{"modTempEff"} fit.
}
\author{ Vito Muggeo}
\seealso{ \code{\link{print.modTempEff}} }
\examples{
  \dontrun{
  #continues from ?tempeff
  anova(o3,o1,test="Chisq") #approximate p-value..
  
  anova(o3,o2,o1,test="Cp")
  
  }
  }
\keyword{ regression}
