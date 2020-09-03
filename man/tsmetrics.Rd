\name{tsmetrics.tsets.estimate}
\alias{tsmetrics}
\alias{tsmetrics.tsets.estimate}
\alias{tsmetrics.tsets.predict}
\title{
Performance Metrics
}
\description{
Performance metrics from an estimated or predicted tsets model.
}
\usage{
\method{tsmetrics}{tsets.estimate}(object, ...)
\method{tsmetrics}{tsets.predict}(object, actual, alpha = 0.1, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate}
}
\item{actual}{
The actual data matched to the dates of the forecasts.
}
\item{alpha}{
The coverage level for distributional forecast metrics.
}
\item{\dots}{
Not currently used.
}
}
\author{
Alexios Galanos
}
