\name{tsprofile}
\alias{tsprofile}
\alias{tsprofile.tsets.estimate}
\title{
Model Simulation Based Profiling
}
\description{
Profiling of model dynamics using simulation/estimation/prediction.
}
\usage{
\method{tsprofile}{tsets.estimate}(object, h = 1, nsim = 100,
seed = NULL, cores = 1, trace = 0, sigma_scale = 1,
solver = "nlminb", autodiff = FALSE, ...)
}
\arguments{
\item{object}{
An object of class tsets.estimate
}
\item{h}{
The forecast horizon on which to evaluate performance metrics.
}
\item{nsim}{
The number of paths to generate.
}
\item{seed}{
An object specifying if and how the random number generator should be initialized.
See the simulate documentation for more details.
}
\item{cores}{
The number of cores to use for parallel estimation of the paths.
}
\item{trace}{
Whether to show a progress bar for the estimation of the paths.
}
\item{sigma_scale}{
A scaling factor for the innovations sd.
}
\item{solver}{
Choice of solver to use for the estimation of the paths.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsetsad package.
}
\item{\dots}{
Not currently used
}
}
\details{
The function profiles an estimated model by simulating and then estimating multiple paths from the
assumed DGP while leaving h values out for prediction evaluation. Each simulated path is equal to
the size of the original dataset plus h additional values, and initialized with the initial state
vector from the model. The resulting output contains the distribution of the MAPE, percent bias (BIAS)
and mean squared log relative error (MSLRE) per horizon h. Since these matrices are of class
\dQuote{tsmodel.distribution} they can be readily plotted with the special purpose \dQuote{plot} function
for this class from the \dQuote{tsmethods} package. Additionally, a data.table matrix is also return with
the distribution of the coefficients from each path estimation.
}
\value{
See details above.
}
\author{
Alexios Galanos
}
