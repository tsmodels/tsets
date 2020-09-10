ets_modelspec <- function(y, model = "AAN", damped = FALSE, power = FALSE, xreg = NULL, frequency = NULL,
                          lambda = NULL, normalized_seasonality = TRUE, fixed_pars = NULL, 
                          scale = FALSE, seasonal_init = "fixed", lambda_lower = 0, 
                          lambda_upper = 1, sampling = NULL, ...)
{
  # 1. Check y
  if  (!is.xts(y)) {
    stop("y must be an xts object")
  }
  if (any(is.na(y))) {
    stop("\nNAs found in y...not allowed.")
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
  if (!is.null(lambda)){
    if(!is.na(lambda) & lambda == 1) lambda <- NULL
  }
  if (substr(model,1,1) == "M") {
    if (!is.null(lambda)) {
      warning("\nMultiplicative error model cannot use a Box Cox transformation (lambda). Setting to NULL.")
      transform <- NULL
    } else{
      transform <- NULL
    }
  } else{
    y_orig <- y
    if (!is.null(lambda)) {
        transform <- box_cox(lambda = lambda, lower = lambda_lower, upper = lambda_upper)
        y <- transform$transform(y = y, frequency = frequency)
        transform$lambda <- attr(y, "lambda")
    } else{
      transform <- NULL
    }
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

  if (scale & model_type(model, power) == 1) {
    scaler <- max(y)
    y <- y/scaler
  } else {
    scaler <- 1
  }
  attr(scaler, "scale") <- scale

  spec$model$model <- validate_model(model)
  spec$model$type <- model_type(model, power)
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
  p <- check_fixed(p, fixed_pars)

  # check fixed parameters
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
  spec$model$parmatrix <- p
  spec$model$power <- power
  spec$model$theta_type_model <- theta_type_model
  spec$model$seasonal_init <- seasonal_init
  class(spec) = c("tsets.spec","tsmodel.spec")
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
