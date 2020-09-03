\name{tsdiagnose}
\alias{tsdiagnose.tsets.estimate}
\alias{tsdiagnose}
\title{
Model Diagnostics
}
\description{
Creates a short summary of model based diagnosticts.
}
\usage{
\method{tsdiagnose}{tsets.estimate}(object, plot = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate}.
}
\item{plot}{
Whether to generate diagnostic plots to accompany summary.
}
\item{\dots}{
Not currently used.
}
}
\value{
A list of tables (printed out and returned insivibly) with
Ljung-Box test for residual autocorrelation, parameter and 
model bounds diagnostics and outlier dates using the Rosner
Test from the \code{EnvStats} package.
}
\author{
Alexios Galanos.
}
