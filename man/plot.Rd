% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.tsets.estimate}
\alias{plot.tsets.estimate}
\alias{plot}
\alias{plot.tsets.simulate}
\alias{plot.tsets.profile}
\title{Object Plots}
\usage{
\method{plot}{tsets.estimate}(x, y = NULL, ...)

\method{plot}{tsets.simulate}(x, y = NULL, ...)

\method{plot}{tsets.profile}(x, y = NULL, type = c("metrics", "coef"), ...)
}
\arguments{
\item{x}{an object of class \dQuote{tsets.estimate}, \dQuote{tsets.simulate}
or \dQuote{tsets.profile}.}

\item{y}{not used.}

\item{...}{additional arguments passed to the underlying plot function.}

\item{type}{type of profile plot for objects of class \dQuote{tsets.profile}.}
}
\description{
Plots for objects generated from the tsets functions.
}
