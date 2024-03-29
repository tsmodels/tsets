#' Model Specification
#'
#' @description Specifies an ETS model prior to estimation.
#' @details The specification object holds the information and data which is
#' then passed to the maximum likelihood estimation routines.
#' @param y an xts vector.
#' @param model the type of model (based on the taxonomy in Hyndman) where the
#' first letter denotes the type of error (E), the second letter the type of
#' trend (T) and the third letter the type of seasonality (S). A value of N
#' denotes none. Models supported are \dQuote{AAA}, \dQuote{AAN}, \dQuote{ANN},
#' \dQuote{ANA}, \dQuote{MMM}, \dQuote{MMN}, \dQuote{MNN}, \dQuote{MNM},
#' \dQuote{MAM} and \dQuote{MAN}. In addition, one may set model equal to
#' \dQuote{Theta}, which is equivalent to a \dQuote{AAN} model with slope
#' parameter (beta) fixed to zero.
#' @param damped whether to include a damped trend.
#' @param power whether to use the power model (only for multiplicative error
#' with additive trend combination), i.e. the \dQuote{MAM} and \dQuote{MAN}
#' models.
#' @param xreg an optional xts matrix of regressors (pre-lagged).
#' @param frequency the frequency of y (if using a seasonal model).
#' @param transformation a valid transformation for y from the \dQuote{tstransform}
#' function in the \dQuote{tsaux} package (currently box-cox or logit are
#' available) applied to additive models only.
#' @param lambda the Box Cox power parameter (lambda). If NA will estimate this
#' using the method of Guerrero.
#' @param lower the lower bound for the transformation.
#' @param upper the upper bound for the transformation.
#' @param normalized_seasonality whether to impose Roberts-McKenzie normalized
#' seasonality.
#' @param fixed_pars a named vector of valid parameter names with values which
#' will be fixed rather than estimated. Valid values are as follows:
#' \itemize{
#' \item \strong{alpha} the adjustment coefficient on the Level component
#' \item \strong{beta} the adjustment coefficient on the Slope component
#' \item \strong{gamma} the adjustment coefficient on the Seasonal component\\
#' \item \strong{phi} the damping parameter
#' \item \strong{theta} the power exponent for the Level and Slope components
#' in the power model
#' \item \strong{delta} the power exponent for the Seasonal component in the power model
#' \item \strong{l0} the initial state value for the Level component
#' \item \strong{b0} the initial state value for the Slope component
#' \item \strong{s0},\ldots, \strong{s[m-1]} e.g. s11, the initial state values for the Seasonal component
#' \item \strong{rho1},\ldots,\strong{rho[k]} e.g. rho12, the coefficients on the regressors
#' }
#' @param scale whether to rescale the data using y/max(y) (only for additive models).
#' This sometimes helps in the optimization.
#' @param seasonal_init whether the initial seasonal states are estimated or
#' fixed (set to a backcast approximation).
#' @param sampling sampling frequency of the dataset. If NULL, will try to
#' identify from the timestamps of y. This is useful for plotting and extending
#' the timestamps in the prediction horizon.
#' @param xreg_init whether to find initial estimates for the regressors with
#' tighter lower and upper bounds. This is only applicable for additive error
#' models.
#' @param ... not used.
#' @return An object of class \dQuote{tsets.spec}.
#' @aliases ets_modelspec
#' @rdname ets_modelspec
#' @export
#'
#'
#'
ets_modelspec <- function(y, model = "AAN", damped = FALSE, power = FALSE, xreg = NULL, frequency = NULL,
                          transformation = "box-cox", lambda = NULL, normalized_seasonality = TRUE,
                          fixed_pars = NULL, scale = FALSE, seasonal_init = "fixed", lower = 0,
                          upper = 1, sampling = NULL, xreg_init = TRUE, ...)
{
  # 1. Check y
  if  (!is.xts(y)) {
    stop("y must be an xts object")
  }
  # 2. Validate model choice
  model <- match.arg(model, choices = c("AAA","AAN","ANN","ANA","MMM","MMN","MNM","MNN","MAM","MAN","Theta"), several.ok = FALSE)

  if (model == "Theta") {
    theta_type_model <- TRUE
    model <- "AAN"
    fixed_pars <- 0
    names(fixed_pars) <- "beta"
  } else{
    theta_type_model <- FALSE
  }
  model_cl <- model_class(model)
  # 3. Validate model inputs
  damped <- as.logical(damped[1])
  power  <- as.logical(power[1])
  normalized_seasonality <- as.logical(normalized_seasonality[1])

  if (substr(model, 1, 1) == "A" & power) {
    warning("\npower model not available with additive error...setting to FALSE")
    power <- FALSE
  }

  if (substr(model,2,2) == "M" & power) {
    warning("\npower model not available with multiplicative trend...setting to FALSE")
    power <- FALSE
  }
  # 4. Check regressors
  xreg <- check_xreg(xreg, index(y))
  # 5. Check transformation
  y_orig <- y
  if (!is.null(transformation)) {
      if (!is.null(lambda) & transformation == "box-cox") {
          if (!is.na(lambda) & lambda == 1) lambda <- NULL
      } else if (transformation != "logit") {
          lambda <- NULL
      }
      if (substr(model,1,1) == "M") {
          if (!is.null(lambda) & !is.null(transformation)) {
              warning("\nMultiplicative error model cannot use a Box Cox or logit transformation. Setting to NULL.")
              transformation <- NULL
              transform <- NULL
              lambda <- NULL
          } else {
            transformation <- NULL
            transform <- NULL
            lambda <- NULL
          }
      } else {
          y_orig <- y
          if (!is.null(lambda) & transformation == "box-cox") {
              if (is.na(lambda)) estimated <- TRUE else estimated <- FALSE
              transform <- tstransform(method = transformation[1], lambda = lambda,
                                       lower = lower, upper = upper)
              y <- transform$transform(y = y, frequency = frequency)
              transform$lambda <- attr(y, "lambda")
              transform$estimated <- estimated
          } else if (transformation == "logit") {
              transform <- tstransform(method = transformation[1], lower = lower, upper = upper)
              y <- transform$transform(y = y)
              transform$estimated <- FALSE
          } else {
              transform <- NULL
          }
      }
  } else {
      transform <- NULL
  }
  if (!is.null(transform)) {
      transform$name <- transformation[1]
      transform$lower <- lower
      transform$upper <- upper
  }
  # 6. Check seasonality
  if (is.null(frequency) & substr(model,3,3) != "N") {
    stop("frequency cannot be NULL when using a seasonal model.")
  }
  spec <- list()
  spec$target$frequency <- frequency

  if (is.null(sampling)) {
    sampling <- sampling_frequency(index(y))
  }
  spec$target$sampling <- sampling

  if (scale & model_type(model, power) != 1) {
    warning("\nscaling can only be used with additive error models. Setting to FALSE.")
    scale <- FALSE
  }

  if (scale & model_type(model, power) == 1 & (!is.null(transform) && transformation == "box-cox")) {
    scaler <- max(y, na.rm = TRUE)
    y <- y/scaler
  } else {
    scaler <- 1
  }
  attr(scaler, "scale") <- scale

  spec$model$model <- validate_model(model)
  spec$model$type <- model_type(model, power)
  spec$model$class <- model_cl
  spec$model$damped <- as.logical(damped)
  spec$model$include_trend <- ifelse(substr(model,2,2) != "N", 1, 0)
  spec$model$include_seasonal <- ifelse(substr(model,3,3) != "N", 1, 0)
  spec$model$include_xreg <- ifelse(!is.null(xreg), 1, 0)
  spec$model$include_damped <- ifelse(spec$model$include_trend == 1 && damped, 1, 0)

  if (power && spec$model$include_seasonal == 1 && normalized_seasonality) {
    warning("\npower model does not support normalized seasonality")
    normalized.seasonality <- FALSE
  }
  spec$model$normalized_seasonality <- ifelse(normalized_seasonality,1,0)

  if (is.null(frequency)) {
    # set to minimum default
    spec$seasonal$frequency <- 4
  } else {
    if (spec$model$include_seasonal == 0) {
      spec$seasonal$frequency <- 4
    } else {
      spec$seasonal$frequency <- frequency
    }
  }

  p <- init_model(model = model, damped = damped, power = power, frequency = spec$seasonal$frequency, y = y, xreg = xreg)
  # check fixed parameters
  p <- check_fixed(p, fixed_pars)
  spec$target$y <- as.numeric(y)
  spec$target$y_orig <- as.numeric(y_orig)
  spec$target$scaler <- scaler
  spec$target$index <- index(y_orig)
  spec$transform <- transform
  spec$model$error <- ifelse(substr(model,1,1) == "A", "Additive", "Multiplicative")
  spec$model$trend <- switch(substr(model,2,2), "A" = "Additive", "M" = "Multiplicative", "N" = "None")
  if (is.null(xreg)) {
    if (spec$model$type == 2) {
      xreg <- matrix(0, ncol = 1, nrow = NROW(y))
    } else{
      xreg <- matrix(0, ncol = 1, nrow = NROW(y))
    }
    include_xreg <- 0
  } else{
    include_xreg <- 1
  }
  spec$xreg$xreg <- coredata(xreg)
  if (!is.null(xreg)) {
    spec$xreg$index <- index(xreg)
  }
  spec$xreg$include_xreg <- include_xreg

  if (theta_type_model) {
    p["beta", c(1,2)] <- 0
  }
  if (spec$model$include_seasonal == 1) {
    if (seasonal_init == "fixed") {
      fixed_s_names <- paste0("s",0:(frequency - 2))
      p[fixed_s_names,"estimate"] <- 0
      p[fixed_s_names,"fixed"] <- 1
    }
  }
  # if alpha, beta or gamma are fixed for models AAA or MAM then we need to impose lower or upper bounds
  if (spec$model$type == 1) {
    if (p["alpha", "fixed"] == 1 & p["beta","estimate"] == 1) {
      p["beta","upper"] <- p["alpha", "init"] - 1e-6
    }
    if (p["alpha", "fixed"] == 1 & p["gamma","estimate"] == 1) {
      p["gamma","upper"] <- 1 - p["alpha", "init"] - 1e-6
    }
    if (p["beta", "fixed"] == 1 & p["alpha","estimate"] == 1) {
      p["alpha","lower"] <- p["beta", "init"] + 1e-6
    }
  }
  if (spec$model$type == 3) {
    if (p["alpha", "fixed"] == 1 & p["beta","estimate"] == 1) {
      p["beta","upper"] <- p["alpha", "init"] - 1e-6
    }
    if (p["beta", "fixed"] == 1 & p["alpha","estimate"] == 1) {
      p["alpha","lower"] <- p["beta", "init"] + 1e-6
    }
  }
  spec$model$parmatrix <- p
  spec$model$power <- power
  spec$model$theta_type_model <- theta_type_model
  spec$model$seasonal_init <- seasonal_init
  class(spec) = c("tsets.spec","tsmodel.spec")
  if (xreg_init) {
    initx <- init_x(spec)
    if (!is.null(initx)) {
      spec$model$parmatrix[grepl("rho",rownames(spec$model$parmatrix)),"init"] <- unname(initx$pars)
      spec$model$parmatrix[grepl("rho",rownames(spec$model$parmatrix)),"lower"] <- unname(initx$lower)
      spec$model$parmatrix[grepl("rho",rownames(spec$model$parmatrix)),"upper"] <- unname(initx$upper)
    }
  }
  return(spec)
}


################################################
# Fixed Parameters
################################################
check_fixed <- function(pars, fixed_pars)
{
  if (is.null(fixed_pars)) {
    return(pars)
  }
  n <- length(fixed_pars)
  fixed_names <- names(fixed_pars)
  # valid values
  pnames <- rownames(pars[which(pars[,"required"] == 1), ])
  for (i in 1:n) {
    if (any(fixed_names[i] %in% pnames)) {
      if (pars[fixed_names[i],"required"] == 1) {
        pars[fixed_names[i],"init"] = fixed_pars[i]
        if (pars[fixed_names[i],"estimate"] == 1) {
          pars[fixed_names[i],"estimate"] <- 0
          pars[fixed_names[i],"fixed"] <- 1
        }
      }
    }
  }
  return(pars)
}


model_class <- function(model)
{
  x <- rbind(
    data.frame(class = 1, models = c("ANN","ANA","AAN","AAA")),
    data.frame(class = 2, models = c("MNN","MNA","MAN","MAA")),
    data.frame(class = 3, models =  c("MNM","MAM")),
    data.frame(class = 4, models = c("MMN","MMM")),
    data.frame(class = 5, models = c("MMA","ANM","AMN","AAM","AMA","AMM")))
  modelc <- x[which(x$models == model),]$class
  return(modelc)
}
