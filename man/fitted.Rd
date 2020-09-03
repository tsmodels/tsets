\name{fitted}
\alias{fitted}
\alias{fitted.tsets.estimate}
\title{
Model Fitted Values
}
\description{
Extract the fitted values from an estimated model.
}
\usage{
\method{fitted}{tsets.estimate}(object, ...)
}
\arguments{
\item{object}{
An object of class tsets.estimate
}
\item{\dots}{
Not currently used
}
}
\value{
An xts vector, back adjusted if a transformation was used (for additive
error models only).
}
\author{
Alexios Galanos
}
