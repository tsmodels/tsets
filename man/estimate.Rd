\name{estimate}
\alias{estimate}
\alias{estimate.tsets.spec}
\title{Model Estimation}
\description{
Estimates a model given a specification object using maximum likelihood.
}
\usage{
\method{estimate}{tsets.spec}(object, solver = "nlminb",
control = list(trace = 0),  autodiff = FALSE, ...)
}
\arguments{
\item{object}{
An object of class tsets.spec.
}
\item{solver}{
One of either \dQuote{solnp}, \dQuote{nlminb} or \dQuote{optim}. The
latter uses the L-BFGS-B algorithm from the lbfgsb3c package. For option
\dQuote{autodiff}, valid solvers are \dQuote{nlminb} and \dQuote{nloptr}.
}
\item{control}{
Solver control parameters.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsetsad package.
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
When autodiff is TRUE with the nloptr solver, the constraints and their jacobian
are explicitly used, whilst a soft barrier constraint is used in the case of
the nlminb solver.
}
\author{
Alexios Galanos
}
