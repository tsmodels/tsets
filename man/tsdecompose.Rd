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
\method{tsdecompose}{tsets.estimate}(object, simplify = FALSE, ...)
\method{tsdecompose}{tsets.predict}(object,  simplify = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate} or \dQuote{tsets.predict}.
}
\item{simplify}{
Whether to return the components as Trend (Level + Slope), Seasonal,
X and Irregular.
}
\item{\dots}{
Not currently used
}
}
\value{
For the estimated object, returns an xts matrix of the state components
(including error). For the predicted object, a list with the
simulated object which includes the distribution matrix of the
state components and their estimated values, inheriting
class \dQuote{tsmodel.predict} (for pretty plotting functionality).
}
\author{
Alexios Galanos
}
