\name{tsmoments}
\alias{tsmoments}
\alias{tsmoments.tsets.estimate}
\title{Analytic Forecast Moments}
\description{
Analytic Forecast Moments
}
\usage{
\method{tsmoments}{tsets.estimate}(object, h = 1, newxreg = NULL,
init_states = NULL, ...)
}
\arguments{
\item{object}{
An object of class tsets.estimate.
}
\item{h}{
The forecast horizon.
}
\item{newxreg}{
The matrix of external regressors in the forecast horizon.
}
\item{init_states}{
An optional vector of states to initialize the forecast. If NULL,
will use the last available state from the estimated model.
}
\item{\dots}{
Not currently used.
}
}
\note{
The function is currently incomplete as it does not provide an option to
backtransform the moments for Box Cox or logit transformed additive models.
}
\author{
Alexios Galanos
}
