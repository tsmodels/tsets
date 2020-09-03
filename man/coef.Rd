\name{coef.tsets.estimate}
\alias{coef}
\alias{coef.tsets.estimate}
\alias{logLik.tsets.estimate}
\alias{logLik}
\alias{AIC.tsets.estimate}
\alias{AIC}
\title{
Estimated Model Extractor Functions
}
\description{
Coefficient, log-likelihood and AIC from an estimated tsets model.
}
\usage{
\method{coef}{tsets.estimate}(object, ...)
\method{logLik}{tsets.estimate}(object, ...)
\method{AIC}{tsets.estimate}(object, ..., k = 2)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate}
}
\item{k}{
(numeric) The penalty per parameter to be used;
the default k = 2 is the classical AIC.
}
\item{\dots}{
Not currently used.
}
}
\author{
Alexios Galanos
}