\name{auto_ets}
\alias{auto_ets}
\title{
Automatic Model Selection
}
\description{
Automatic model selection based on choice of criteria.
}
\usage{
auto_ets(y, xreg = NULL, transformation = NULL, lambda = NULL,
    lower = 0, upper = 1, metric = "AIC", frequency = NULL, normalized_seasonality = TRUE,
    additive_only = FALSE, cores = NULL, solver = "nlminb", control = list(trace = 0,
        maxit = 1000), power_model = FALSE, include_damped = TRUE,
    verbose = FALSE, retain = 1, scale = FALSE, seasonal_init = "fixed",
    autodiff = FALSE, ...)
}
\arguments{
\item{y}{
An xts vector.
}
\item{xreg}{
An optional xts matrix of regressors
}
\item{transformation}{
A valid transformation for y from the \dQuote{tstransform} function in the
\dQuote{tsaux} package (current box-cox or logit are available) applied to
additive models only.
}
\item{lambda}{
The Box Cox power parameter (lambda). If NA will estimate this
using the method of Guerrero.
}
\item{lower}{
The lower bound for the transformation.
}
\item{upper}{
The upper bound for the transformation.
}
\item{metric}{
The selection metric to use. Valid metrics are \sQuote{AIC}, \sQuote{BIC}, \sQuote{AICc},
\sQuote{MASE} and \sQuote{MAPE}. If lambda is not NULL, then those models which
admit a Box-Cox transformation (additive error models) will not be comparable with the
other models.
}
\item{frequency}{
The frequency of y (if using a seasonal model).
}
\item{normalized_seasonality}{
Whether to impose Roberts-McKenzie normalized seasonality.
}
\item{additive_only}{
Whether to limit to additive models only.
}
\item{cores}{
The number of cores to use for parallel evaluation.
}
\item{solver}{
The solver to use for estimation.
}
\item{control}{
The solver control parameters.
}
\item{power_model}{
Wether to include the power MAM models.
}
\item{include_damped}{
Whether to include damped models in the selection.
}
\item{verbose}{
Whether to generate verbose output.
}
\item{retain}{
The models to retain (ranked by metric).
}
\item{scale}{
Optional scaling for additive models.
}
\item{seasonal_init}{
Whether the initial seasonal states are estimated or
fixed (set to a backcast approximation).
}
\item{autodiff}{
Whether to use automatic differentiation (see \code{\link{estimate.tsets.spec}}).
}
\item{\dots}{
Not currently used.
}
}
\details{
The function estimates a number of different models from the class
of valid ets models in the packages, with and without dampening,
with and without a power term (for MAM and MAN models).
}
\value{
An object of class \dQuote{tsets.estimate} which also inherits
class \dQuote{tsets.select}
}
\author{
Alexios Galanos
}
