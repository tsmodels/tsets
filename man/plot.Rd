\name{plot}
\alias{plot}
\alias{plot.tsets.estimate}
\alias{plot.tsets.simulate}
\alias{plot.tsets.profile}
\title{
Model Plots
}
\description{
Plots for objects generated from the tsets functions.
}
\usage{
\method{plot}{tsets.estimate}(x, y = NULL, ...)
\method{plot}{tsets.simulate}(x, y = NULL, ...)
\method{plot}{tsets.profile}(x, y = NULL, type = c("metrics", "coef"), ...)
}
\arguments{
\item{x}{
An object of class \dQuote{tsets.estimate}.
}
\item{y}{
Not used.
}
\item{type}{
The type of profile plot.
}
\item{\dots}{
Additional arguments passed to the underlying plot function.
}
}
\author{
Alexios Galanos
}
