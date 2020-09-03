\name{estimate}
\alias{estimate}
\alias{estimate.tsets.spec}
\title{Model Estimation}
\description{
Estimates a model given a specification object using maximum likelihood.
}
\usage{
\method{estimate}{tsets.spec}(object, solver = "nlminb",
control = list(trace = 0),  ...)
}
\arguments{
\item{object}{
An object of class tsets.spec.
}
\item{solver}{
One of either \dQuote{solnp}, \dQuote{nlminb} or \dQuote{optim}. The
latter uses the L-BFGS-B algorithm from the lbfgsb3c package.
}
\item{control}{
Solver control parameters.
}
\item{\dots}{
Not currently used.
}
}
\details{
The maximum likelihood estimation uses bound constraints for some of the parameters
as described in the vignette of the package. Additionally, for parameters which
are constrained to be less than another parameter (e.g. beta<alpha), a simple
barrier approach is adopted which adjust the previous likelihood value upwards by
some fixed percentage of that value during the minimization.
The observation variance is not directly estimated but instead concentrated out
of the likelihood.
}
\author{
Alexios Galanos
}
