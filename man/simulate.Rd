\name{simulate.tsets.estimate}
\alias{simulate.tsets.estimate}
\alias{simulate}
\title{
Model Simulation
}
\description{
Simulates paths from an ETS model with optional user specified values.
}
\usage{
\method{simulate}{tsets.estimate}(object, nsim = 1, seed = NULL, 
h = NULL, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, 
innov = NULL, sigma_scale = 1, pars = coef(object), ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate}.
}
\item{nsim}{
The number of paths per complete set of time steps (h).
}
\item{seed}{
An object specifying if and how the random number generator should be
initialized (\sQuote{seeded}). Either NULL or an integer that will be used
in a call to set.seed before simulating the response vectors. If set,
the value is saved as the "seed" attribute of the returned value.
The default, NULL will not change the random generator state, and
return .Random.seed as the \dQuote{seed} attribute in the returned object.
}
\item{h}{
The number of time steps to simulate paths for. If this is NULL, it will
use the same number of periods as in the original series.
}
\item{newxreg}{
An optional matrix of regressors to use for the simulation if xreg
was used in the estimation. If NULL and the estimated object had
regressors, and h was also set to NULL, then the original regressors
will be used.
}
\item{sim_dates}{
An optional vector of simulation dates equal to h. If NULL will
use the implied periodicity of the data to generate a regular
sequence of dates after the first available date in the data.
}
\item{bootstrap}{
Whether to bootstrap the innovations from the estimated object
by resampling from the empirical distribution.
}
\item{innov}{
This is an optional vector of uniform innovations which will be
translated to regular innovations using the appropriate
distribution quantile function and model standard deviation.
The length of this vector should be equal to nsim x horizon.
}
\item{sigma_scale}{
An optional scalar which will scale the standard deviation of
the innovations (useful for profiling under different assumptions).
This applied to all approaches (parametric, bootstrap and custom
innov supplied values).
}
\item{pars}{
An optional named vector of model coefficients which ovveride the
estimated coefficients. No checking is currently performed on the
adequacy of these coefficients.
}
\item{\dots}{
Not currently used.
}
}
\details{
It is required that an initial object of class \dQuote{tsets.estimate}
be passed to this function, but the parameters and initial states can
be overriden by passing them as named values in the \sQuote{pars} 
argument. The default is to initialize the states from the seed
states used in the estimated object, with h equal to the length
of the original series (default for NULL h).\cr
Innovations for the simulation can either be parametric (Normal
for additive or truncated Normal for multiplicative error models),
based on the estimated residuals (bootstrap argument) or use 
supplied set of uniform random number (which are then translated
into Normal or truncated Normal using sd equal to the model sigma
and optionally scaled by sigma_scale).
}
\value{
An object of class \dQuote{tsets.simulate} with slots for the
simulated series and states.
}
\author{
Alexios Galanos
}
