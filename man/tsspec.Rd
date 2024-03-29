% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{tsspec.tsets.estimate}
\alias{tsspec.tsets.estimate}
\alias{tsspec}
\title{Model Specification Extractor}
\usage{
\method{tsspec}{tsets.estimate}(object, y = NULL, lambda = NULL, xreg = NULL, ...)
}
\arguments{
\item{object}{an object of class \dQuote{tsets.estimate}.}

\item{y}{an optional new xts vector.}

\item{lambda}{an optional lambda parameter for the Box Cox transformation (if
previously used).}

\item{xreg}{an optional matrix of regressors.}

\item{...}{not currently used.}
}
\value{
An object of class \dQuote{tsets.spec}.
}
\description{
Extracts a model specification (class \dQuote{tsets.spec}) from
an object of class \dQuote{tsets.estimate}.
}
\note{
This function is used by other functions in the package such as the
backtest which requires rebuilding the specification for each re-estimation
step with updated data but keeping all else equal.
}
