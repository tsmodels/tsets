% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specification.R
\name{ets_modelspec}
\alias{ets_modelspec}
\title{Model Specification}
\usage{
ets_modelspec(
  y,
  model = "AAN",
  damped = FALSE,
  power = FALSE,
  xreg = NULL,
  frequency = NULL,
  transformation = "box-cox",
  lambda = NULL,
  normalized_seasonality = TRUE,
  fixed_pars = NULL,
  scale = FALSE,
  seasonal_init = "fixed",
  lower = 0,
  upper = 1,
  sampling = NULL,
  xreg_init = TRUE,
  ...
)
}
\arguments{
\item{y}{an xts vector.}

\item{model}{the type of model (based on the taxonomy in Hyndman) where the
first letter denotes the type of error (E), the second letter the type of
trend (T) and the third letter the type of seasonality (S). A value of N
denotes none. Models supported are \dQuote{AAA}, \dQuote{AAN}, \dQuote{ANN},
\dQuote{ANA}, \dQuote{MMM}, \dQuote{MMN}, \dQuote{MNN}, \dQuote{MNM},
\dQuote{MAM} and \dQuote{MAN}. In addition, one may set model equal to
\dQuote{Theta}, which is equivalent to a \dQuote{AAN} model with slope
parameter (beta) fixed to zero.}

\item{damped}{whether to include a damped trend.}

\item{power}{whether to use the power model (only for multiplicative error
with additive trend combination), i.e. the \dQuote{MAM} and \dQuote{MAN}
models.}

\item{xreg}{an optional xts matrix of regressors (pre-lagged).}

\item{frequency}{the frequency of y (if using a seasonal model).}

\item{transformation}{a valid transformation for y from the \dQuote{tstransform}
function in the \dQuote{tsaux} package (currently box-cox or logit are
available) applied to additive models only.}

\item{lambda}{the Box Cox power parameter (lambda). If NA will estimate this
using the method of Guerrero.}

\item{normalized_seasonality}{whether to impose Roberts-McKenzie normalized
seasonality.}

\item{fixed_pars}{a named vector of valid parameter names with values which
will be fixed rather than estimated. Valid values are as follows:
\itemize{
\item \strong{alpha} the adjustment coefficient on the Level component
\item \strong{beta} the adjustment coefficient on the Slope component
\item \strong{gamma} the adjustment coefficient on the Seasonal component\\
\item \strong{phi} the damping parameter
\item \strong{theta} the power exponent for the Level and Slope components
in the power model
\item \strong{delta} the power exponent for the Seasonal component in the power model
\item \strong{l0} the initial state value for the Level component
\item \strong{b0} the initial state value for the Slope component
\item \strong{s0},\ldots, \strong{s[m-1]} e.g. s11, the initial state values for the Seasonal component
\item \strong{rho1},\ldots,\strong{rho[k]} e.g. rho12, the coefficients on the regressors
}}

\item{scale}{whether to rescale the data using y/max(y) (only for additive models).
This sometimes helps in the optimization.}

\item{seasonal_init}{whether the initial seasonal states are estimated or
fixed (set to a backcast approximation).}

\item{lower}{the lower bound for the transformation.}

\item{upper}{the upper bound for the transformation.}

\item{sampling}{sampling frequency of the dataset. If NULL, will try to
identify from the timestamps of y. This is useful for plotting and extending
the timestamps in the prediction horizon.}

\item{xreg_init}{whether to find initial estimates for the regressors with
tighter lower and upper bounds. This is only applicable for additive error
models.}

\item{...}{not used.}
}
\value{
An object of class \dQuote{tsets.spec}.
}
\description{
Specifies an ETS model prior to estimation.
}
\details{
The specification object holds the information and data which is
then passed to the maximum likelihood estimation routines.
}
