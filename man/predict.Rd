\name{predict}
\alias{predict}
\alias{predict.tsets.estimate}
\title{
Model Prediction
}
\description{
Prediction function for class \dQuote{tsets.estimate}.
}
\usage{
\method{predict}{tsets.estimate}(object, h = 12, newxreg = NULL,
nsim = 1000, drop_na = TRUE, drop_negative = FALSE, redraw = TRUE,
forc_dates = NULL, innov = NULL, custom_slope = NULL, ...)
}
\arguments{
\item{object}{
An object of class tsets.estimate
}
\item{h}{
The forecast horizon.
}
\item{newxreg}{
The matrix of external regressors in the forecast horizon.
}
\item{nsim}{
The number of simulations to use for generating the simulated
predictive distribution.
}
\item{drop_na}{
Drop any simulations that result in NA or NaN values.
}
\item{drop_negative}{
Drop any simulations that result in negative values.
}
\item{redraw}{
Whether to re-draw from the distribution if any NA, NaN or
negative values resulted in a matrix with less than nsim
simulations.
}
\item{forc_dates}{
An optional vector of forecast dates equal to h. If NULL will
use the implied periodicity of the data to generate a regular
sequence of dates after the last available date in the data.
}
\item{innov}{
This is an optional vector of uniform innovations which will be 
translated to regular innovations using the appropriate 
distribution quantile function and model standard deviation. 
The length of this vector should be equal to nsim x horizon.
}
\item{custom_slope}{
Either a vector of length equal to the horizon or 1. This will be
used to ovveride the slope state with a user provided set of 
values (or strong views on growth). Only allowed for AA or MM
type models (i.e. no MA type).
}
\item{\dots}{
Not currently used
}
}
\details{
Like all models in the ts framework, prediction is done by simulating
h-steps ahead in order to build a predictive distribution.
}
\value{
An object of class \dQuote{tsets.predict} which also inherits 
\dQuote{tsmodel.predict}, with slots for the simulated prediction distribution, 
the original series (as a zoo object), the original specification object and 
the mean forecast. The predictive distribution is back transformed if lambda was
not set to NULL in the specification.
}
\author{
Alexios Galanos
}
