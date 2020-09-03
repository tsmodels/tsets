\name{tsdecompose}
\alias{tsdecompose}
\alias{tsdecompose.tsets.estimate}
\alias{tsdecompose.tsets.predict}
\title{
Model Decomposition
}
\description{
Decomposes the estimated model or prediction into its component parts.
}
\usage{
\method{tsdecompose}{tsets.estimate}(object, ...)
\method{tsdecompose}{tsets.predict}(object, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate} or \dQuote{tsets.predict}.
}
\item{\dots}{
Not currently used
}
}
\value{
For the estimated object, returns a list of xts vectors with the fitted and state components (including error).
For the predicted object, a list with the simulated distribution matrix of the optional forecasted and state
components which also inherits class \dQuote{tsmodel.distribution} (for pretty plotting functionality).
}
\author{
Alexios Galanos
}
