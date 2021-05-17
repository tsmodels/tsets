\name{summary}
\alias{summary}
\alias{summary.tsets.estimate}
\title{
Model Estimation Summary
}
\description{
Summary method for class \dQuote{tsets.estimate}
}
\usage{
\method{summary}{tsets.estimate}(object, digits=4, ...)
}
\arguments{
\item{object}{
An object of class tsets.estimate
}
\item{digits}{
Integer, used for number formatting. Optionally, to avoid scientific notation,
set \sQuote{options(scipen=999)}.
}
\item{\dots}{
Not currently used
}
}
\value{
A printout of the parameter summary, model type and some model metrics. When
estimated using autodiff, the standard errors, t-values and p-values
will also be printed. In this case, if the parameters are close to their
upper or lower bounds then it is very likely that these values will be NaN.
}
\author{
Alexios Galanos
}
