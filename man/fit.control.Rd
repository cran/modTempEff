\name{fit.control}
\alias{fit.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary Function for Controlling Model Fitting }
\description{
       Auxiliary function as user interface for fitting. Typically
     only used when calling 'tempeff' or 'tempeff.fit'.
}
\usage{
fit.control(toll = 1e-04, display = FALSE, it.max = 20, 
  GLM = FALSE, maxit.inner = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toll}{ positive convergence tolerance. }
  \item{display}{ logical indicating if deviance should be printed for each iteration. }
  \item{it.max}{ integer giving the maximal number of iterations. }
  \item{GLM}{ logical indicating if at each iteration a GLM (rather than a GAM) has to be fitted. }
  \item{maxit.inner}{ integer giving the maximal number of inner iterations. }
}
\details{
  Fitting of Constrained Distributed Lag Model is performed via iterative estimate of proper
  Generalized Additive (or Linear when \code{GLM=TRUE}) models. The algorithm stops when 
  the relative increase in deviance is smaller than \code{toll} or when the maximum number of 
  iterations \code{it.max} is attained. The maximum number of (inner) iterations to estimate the working 
  GLM at each (outer) iteration is controlled by \code{maxit.inner}. 
}
\value{
  A list with the arguments as components to be used by 'tempeff' or 'tempeff.fit'.
}
\author{ Vito Muggeo }

\seealso{\code{\link{tempeff}}}

\keyword{ regression }
