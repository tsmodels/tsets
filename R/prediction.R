predict.tsets.estimate <- function(object, h = 12, newxreg = NULL, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL,
                                   init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
{
  if (!is.null(forc_dates)) {
    if (h != length(forc_dates)) stop("\nforc_dates must have length equal to h")
  }
  if (!is.null(custom_slope)) {
    if (substr(object$model$setup$model,2,2) == "N") {
        custom_slope <- NULL
        warnings("\ncustom_slope passed to a model without a slope component")
    }
    if (substr(object$model$setup$model,1,2) == "MA") {
      custom_slope <- NULL
      stop("\ncustom_slope passed to an MA type model. Only MM or AA type admit custom_slope.")
    }
  }
  if (is.null(innov_type)) innov_type <- "q"
  innov_type <- match.arg(innov_type, c("q","z"))
  if (object$spec$xreg$include_xreg == 0) {
    newxreg <- NULL
    if (is.null(forc_dates)) {
      forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
    }
  } else {
    if (!is.null(newxreg)) {
      forc_dates <- index(newxreg)
    } else {
      if (is.null(forc_dates)) {
        forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
      }
      warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
      newxreg <- xts(matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h), forc_dates)
      colnames(newxreg) <- colnames(object$spec$xreg$xreg)
    }
  }
  if (!is.null(sigma_scale)) {
    sigma_scale <- as.numeric(sigma_scale)
    if (any(sigma_scale <= 0)) stop("\nsigma_scale must be strictly positive")
    if (length(sigma_scale) == 1) sigma_scale <- rep(sigma_scale, h)
    if (length(sigma_scale) != h) stop("\nsigma_scale must be of length h or 1 (recycled to h)")
  }
  if (!is.null(init_states)) {
    if (length(as.vector(init_states)) != NCOL(object$model$states)) {
      stop(paste0("\ninit_states must be a vector of length ", NCOL(object$model$states)))
    } else {
      init_states <- matrix(as.numeric(init_states), nrow = 1, ncol = NCOL(object$model$states))
    }
  }
  if (!is.null(innov)) {
    requiredn <- h * nsim
    if (length(innov) != requiredn) {
      stop("\nlength of innov must be nsim x h")
    }
    # check that the innovations are uniform samples (from a copula)
    if (innov_type == "q") {
      if (any(innov < 0 | innov > 1 )) {
        stop("\ninnov must be >0 and <1 (uniform samples) for innov_type = 'q'")
      }
      if (any(innov == 0)) innov[which(innov == 0)] <- 1e-12
      if (any(innov == 1)) innov[which(innov == 1)] <- (1 - 1e-12)
    }
    innov <- matrix(innov, h, nsim)
  }
  # run simulation
  zList <- forecast_simulation(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov = innov,
                               custom_slope = custom_slope, init_states = init_states,
                               innov_type = innov_type, sigma_scale = sigma_scale, ...)
  fcast_dist <- zList$distribution
  tmp <- forecast_backtransform(fcast_dist, object$spec$transform)
  fcast_mean <- zoo(tmp$mean, forc_dates)
  fcast_dist <- tmp$dist
  colnames(fcast_dist) <- colnames(zList$distribution)
  # collect and return
  zList$mean <- fcast_mean
  zList$distribution <- fcast_dist
  class(zList$distribution) <- c("tsets.distribution","tsmodel.distribution")
  attr(zList$distribution, "date_class") <- attr(object$spec$target$sampling, "date_class")
  class(zList) <- c("tsets.predict","tsmodel.predict")
  return(zList)
}


forecast_aaa_cpp <- function(object, newxreg = NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
{
  if (is.null(forc_dates)) {
    forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
  }
  coefficient <- object$model$setup$parmatrix[,1]
  frequency <- object$spec$seasonal$frequency
  if (!is.null(custom_slope) & object$model$setup$include_trend == 1) {
    custom_flag <- 1
  } else {
    custom_flag <- 0
  }
  model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim, custom_flag)

  # last state
  if (!is.null(init_states)) {
    stateinit <- init_states
    colnames(stateinit) <- colnames(object$model$states)
  } else {
    stateinit <- tail(object$model$states, 1)
  }
  pars <- rep(0, 6)
  pars[1] <- stateinit[1,"Level"]
  if (model[1] == 1) pars[2] <- stateinit[1,"Trend"]
  pars[3] <- coefficient["alpha"]
  pars[4] <- coefficient["beta"]
  pars[5] <- coefficient["gamma"]
  pars[6] <- coefficient["phi"]
  pars <- unname(as.numeric(pars))

  if (model[2] == 1) {
    s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
  } else {
    s0 <- rep(0, frequency)
  }
  s0 <- unname(as.numeric(s0))

  if (object$model$setup$include_xreg == 1) {
    k <- ncol(object$spec$xreg$xreg)
    if (ncol(newxreg) != k) {
      stop("\nNCOL newxreg not equal to number of columns of fitted xreg!")
    }
    if (NROW(newxreg) != h) {
      stop("\nNROW newxreg not equal to forecast horizon h")
    }
    rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
    xregf <- coredata(newxreg) %*% rho
  } else {
    xregf <- rep(0, h)
  }
  xregf <- c(0, as.numeric(xregf))
  if (!is.null(innov)) {
    if (innov_type == "q") {
      E <- t(qnorm(innov, mean = 0, sd = coefficient["sigma"]))
    } else {
      E <- matrix(innov * coefficient["sigma"], nrow = nsim, ncol = h)
    }
    E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
  } else {
    E <- matrix( rnorm(nsim*(h + 1), 0, coefficient["sigma"]), nsim, h + 1 )
  }
  if (!is.null(custom_slope) & model[1] == 1) {
      if (length(custom_slope) == h) {
        B = matrix(c(pars[2], custom_slope), ncol = h + 1, nrow = nsim, byrow = TRUE)
      } else if (length(custom_slope) == 1) {
        B = matrix(c(pars[2], rep(custom_slope, h)), ncol = h + 1, nrow = nsim, byrow = TRUE)
      } else {
        stop("\ncustom_slope: provide a vector of length equal to either h or 1")
      }
  } else {
    B <- matrix(0, ncol = h + 1, nrow = nsim)
  }
  out <- simulate_aaa(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xregf, slope_overide_ = B)

  Y <- out$Simulated[,-1,drop = FALSE]
  if (!is.null(sigma_scale)) {
    mu <- colMeans(Y)
    Y <- sweep(Y, 2, mu, "-")
    Y <- sweep(Y, 2, sigma_scale, "*")
    Y <- sweep(Y, 2, mu, "+")
  }
  Level <- out$Level
  Level <- Level[,1:(ncol(Level) - 1), drop = FALSE]
  Slope <- out$Slope
  Slope <- Slope[,1:(ncol(Slope) - 1), drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,1:(ncol(Seasonal) - 1), drop = FALSE]
  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)
  return(zList)
}

forecast_mmm_cpp <- function(object, newxreg = NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
{
  if (is.null(forc_dates)) {
    forc_dates <- future_dates(tail(object$spec$target$index,1),frequency = object$spec$target$frequency, n = h)
  }
  coefficient <- object$model$setup$parmatrix[,1]
  frequency <- object$spec$seasonal$frequency
  frequency <- object$spec$seasonal$frequency
  if (!is.null(custom_slope) & object$model$setup$include_trend == 1) {
    custom_flag <- 1
  } else {
    custom_flag <- 0
  }
  model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim, custom_flag)
  # last state
  if (!is.null(init_states)) {
    stateinit <- init_states
    colnames(stateinit) <- colnames(object$model$states)
  } else {
    stateinit <- tail(object$model$states, 1)
  }
  pars <- rep(0, 6)
  pars[1] <- stateinit[1,"Level"]
  if (model[1] == 1) pars[2] <- stateinit[1,"Trend"] else pars[2] <- 1
  pars[3] <- coefficient["alpha"]
  pars[4] <- coefficient["beta"]
  pars[5] <- coefficient["gamma"]
  pars[6] <- coefficient["phi"]
  pars <- unname(as.numeric(pars))

  if (model[2] == 1) {
    s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
  } else{
    s0 <- rep(1, frequency)
  }
  s0 <- unname(as.numeric(s0))

  if (object$model$setup$include_xreg == 1) {
    k <- ncol(object$spec$xreg$xreg)

    if (ncol(newxreg) != k) {
      stop("\nNCOL newxreg not equal to number of columns of fitted xreg!")
    }
    if (NROW(newxreg) != h) {
      stop("\nNROW newxreg not equal to forecast horizon h")
    }

    rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
    xregf <- coredata(newxreg) %*% rho
  } else {
    xregf <- rep(0, h)
  }
  xregf <- c(0, as.numeric(xregf))

  if (!is.null(innov)) {
    if (innov_type == "q") {
      E <- t(matrix(qtruncnorm(innov, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim))
      E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
    } else {
      E <- pnorm(innov, 0, 1)
      E <- matrix(qtruncnorm(E, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim)
    }
    E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
  } else {
    E <- matrix(tsaux:::rtruncnorm(nsim*(h + 1), mu = 0, sigma = coefficient["sigma"], lb = -1), nsim, h + 1)
  }
  if (!is.null(custom_slope) & model[1] == 1) {
    if (length(custom_slope) == h) {
      B = matrix(c(pars[2], custom_slope), ncol = h + 1, nrow = nsim, byrow = TRUE)
    } else if (length(custom_slope) == 1) {
      B = matrix(c(pars[2], rep(custom_slope, h)), ncol = h + 1, nrow = nsim, byrow = TRUE)
    } else {
      stop("\ncustom_slope: provide a vector of length equal to either h or 1")
    }
  } else {
    B <- matrix(0, ncol = h + 1, nrow = nsim)
  }
  out <- simulate_mmm(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xregf, slope_overide_ = B)
  Y <- out$Simulated[,-1,drop = FALSE]
  if (!is.null(sigma_scale)) {
    mu <- colMeans(Y)
    Y <- sweep(Y, 2, mu, "-")
    Y <- sweep(Y, 2, sigma_scale, "*")
    Y <- sweep(Y, 2, mu, "+")
  }
  Level <- out$Level
  Level <- Level[,1:(ncol(Level) - 1), drop = FALSE]
  Slope <- out$Slope
  Slope <- Slope[,1:(ncol(Slope) - 1), drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,1:(ncol(Seasonal) - 1), drop = FALSE]
  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)
  return(zList)
}


forecast_mam_cpp <- function(object, newxreg=NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
{
  if (is.null(forc_dates)) {
    forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$frequency, n = h)
  }
  coefficient <- object$model$setup$parmatrix[,1]
  frequency <- object$spec$seasonal$frequency
  if (!is.null(custom_slope) & object$model$setup$include_trend == 1) {
    custom_flag <- 1
  } else {
    custom_flag <- 0
  }
  model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim, custom_flag)

  # last state
  if (!is.null(init_states)) {
    stateinit <- init_states
    colnames(stateinit) <- colnames(object$model$states)
  } else {
    stateinit <- tail(object$model$states, 1)
  }
  pars <- rep(0, 6)
  pars[1] <- stateinit[1,"Level"]
  if (model[1] == 1) pars[2] <- stateinit[1,"Trend"]
  pars[3] <- coefficient["alpha"]
  pars[4] <- coefficient["beta"]
  pars[5] <- coefficient["gamma"]
  pars[6] <- coefficient["phi"]
  pars <- unname(as.numeric(pars))

  if (model[2] == 1) {
    s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
  } else {
    s0 <- rep(1, frequency)
  }
  s0 <- unname(as.numeric(s0))

  if (object$model$setup$include_xreg == 1) {
    k <- ncol(object$spec$xreg$xreg)

    if (ncol(newxreg) != k) {
      stop("\nNCOL newxreg not equal to number of columns of fitted xreg!")
    }
    if (NROW(newxreg) != h) {
      stop("\nNROW newxreg not equal to forecast horizon h")
    }

    rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
    xregf <- coredata(newxreg) %*% rho
  } else {
    xregf <- rep(0, h)
  }
  xregf <- c(0, as.numeric(xregf))
  if (!is.null(innov)) {
    if (innov_type == "q") {
      E <- t(matrix(qtruncnorm(innov, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim))
      E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
    } else {
      E <- pnorm(innov, 0, 1)
      E <- matrix(qtruncnorm(E, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim)
    }
    E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
  } else {
    E <- matrix(tsaux:::rtruncnorm(nsim*(h + 1), mu = 0, sigma = coefficient["sigma"], lb = -1), nsim, h + 1)
  }

  if (!is.null(custom_slope) & model[1] == 1) {
    if (length(custom_slope) == h) {
      B = matrix(c(pars[2], custom_slope), ncol = h + 1, nrow = nsim, byrow = TRUE)
    } else if (length(custom_slope) == 1) {
      B = matrix(c(pars[2], rep(custom_slope, h)), ncol = h + 1, nrow = nsim, byrow = TRUE)
    } else {
      stop("\ncustom_slope: provide a vector of length equal to either h or 1")
    }
  } else {
    B <- matrix(0, ncol = h + 1, nrow = nsim)
  }
  out <- simulate_mam(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xregf, slope_overide_ = B)
  Y <- out$Simulated[,-1,drop = FALSE]
  if (!is.null(sigma_scale)) {
    mu <- colMeans(Y)
    Y <- sweep(Y, 2, mu, "-")
    Y <- sweep(Y, 2, sigma_scale, "*")
    Y <- sweep(Y, 2, mu, "+")
  }
  Level <- out$Level
  Level <- Level[,1:(ncol(Level) - 1), drop = FALSE]
  Slope <- out$Slope
  Slope <- Slope[,1:(ncol(Slope) - 1), drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,1:(ncol(Seasonal) - 1), drop = FALSE]

  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)

  return(zList)
}

forecast_powermam_cpp <- function(object, newxreg = NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
{
  if (is.null(forc_dates)) {
    forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$frequency, n = h)
  }
  coefficient <- object$model$setup$parmatrix[,1]
  frequency <- object$spec$seasonal$frequency
  if (!is.null(custom_slope) & object$model$setup$include_trend == 1) {
    custom_flag <- 1
  } else {
    custom_flag <- 0
  }

  model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim, custom_flag)

  # last state
  if (!is.null(init_states)) {
    stateinit <- init_states
    colnames(stateinit) <- colnames(object$model$states)
  } else {
    stateinit <- tail(object$model$states, 1)
  }
  pars <- rep(0, 8)
  pars[1] <- stateinit[1,"Level"]
  if (model[1] == 1) pars[2] <- stateinit[1,"Trend"]
  pars[3] <- coefficient["alpha"]
  pars[4] <- coefficient["beta"]
  pars[5] <- coefficient["gamma"]
  pars[6] <- coefficient["phi"]
  pars[7] <- coefficient["theta"]
  pars[8] <- coefficient["delta"]

  pars <- unname(as.numeric(pars))

  if (model[2] == 1) {
    s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
  } else {
    s0 <- rep(1, frequency)
  }
  s0 <- unname(as.numeric(s0))

  if (object$model$setup$include_xreg == 1) {
    k <- ncol(object$spec$xreg$xreg)
    if (ncol(newxreg) != k) {
      stop("\nNCOL newxreg not equal to number of columns of fitted xreg!")
    }
    if (NROW(newxreg) != h) {
      stop("\nNROW newxreg not equal to forecast horizon h")
    }

    rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
    xregf <- coredata(newxreg) %*% rho
  } else {
    xregf <- rep(0, h)
  }

  xregf <- c(0, as.numeric(xregf))
  if (!is.null(innov)) {
    if (innov_type == "q") {
      E <- t(matrix(qtruncnorm(innov, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim))
      E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
    } else {
      E <- pnorm(innov, 0, 1)
      E <- matrix(qtruncnorm(E, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim)
    }
    E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
  } else {
    E <- matrix(tsaux:::rtruncnorm(nsim*(h + 1), mu = 0, sigma = coefficient["sigma"], lb = -1), nsim, h + 1)
  }
  if (!is.null(custom_slope) & model[1] == 1) {
    if (length(custom_slope) == h) {
      B = matrix(c(pars[2], custom_slope), ncol = h + 1, nrow = nsim, byrow = TRUE)
    } else if (length(custom_slope) == 1) {
      B = matrix(c(pars[2], rep(custom_slope, h)), ncol = h + 1, nrow = nsim, byrow = TRUE)
    } else {
      stop("\ncustom_slope: provide a vector of length equal to either h or 1")
    }
  } else {
    B <- matrix(0, ncol = h + 1, nrow = nsim)
  }
  #
  out <- simulate_powermam(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xregf, slope_overide_ = B)
  Y <- out$Simulated[,-1,drop = FALSE]
  if (!is.null(sigma_scale)) {
    mu <- colMeans(Y)
    Y <- sweep(Y, 2, mu, "-")
    Y <- sweep(Y, 2, sigma_scale, "*")
    Y <- sweep(Y, 2, mu, "+")
  }
  Level <- out$Level
  Level <- Level[,1:(ncol(Level) - 1), drop = FALSE]
  Slope <- out$Slope
  Slope <- Slope[,1:(ncol(Slope) - 1), drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,1:(ncol(Seasonal) - 1), drop = FALSE]
  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)
  return(zList)
}


wrap_forecast_output <- function(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)
{
  colnames(Y)         <- as.character(forc_dates)
  colnames(Level)     <- as.character(forc_dates)
  colnames(Slope)     <- as.character(forc_dates)
  colnames(Seasonal)  <- as.character(forc_dates)
  names(xregf)        <- as.character(forc_dates)
  colnames(E)        <- as.character(forc_dates)
  date_class <- attr(object$spec$target$sampling, "date_class")
  dist_classes <- c("tsets.distribution", "tsmodel.distribution")
  class(Level) <- dist_classes
  attr(Level, "date_class") <- date_class
  class(Y)     <- dist_classes
  attr(Y, "date_class") <- date_class
  class(E)     <- dist_classes
  attr(E, "date_class") <- date_class
  tsx <- tsdecompose(object)
  Level <- list(distribution = Level, original_series = tsx$Level)
  class(Level) <- "tsmodel.predict"
  if (model[1] == 0) {
    Slope <- NULL
  } else {
    class(Slope) <- dist_classes
    attr(Slope, "date_class") <- date_class
    Slope <- list(distribution = Slope, original_series = tsx$Slope)
    class(Slope) <- "tsmodel.predict"
  }
  if (model[2] == 0) {
    Seasonal <- NULL
  } else{
    class(Seasonal) <- dist_classes
    attr(Seasonal, "date_class") <- date_class
    Seasonal <- list(distribution = Seasonal, original_series = tsx$Seasonal)
    class(Seasonal) <- "tsmodel.predict"
  }
  if (object$model$setup$include_xreg == 0) {
    xregf <- NULL
  }
  E <- list(distribution = E, original_series = residuals(object, raw = TRUE))
  class(E) <- "tsmodel.predict"
  decomposition <- list(Level = Level, Slope = Slope, Seasonal = Seasonal, X = xregf, Error = E)
  zList <- list(distribution = Y, original_series = zoo(object$spec$target$y_orig, object$spec$target$index), h = model[3],
                spec = object$spec, decomposition = decomposition)
  class(zList) <- dist_classes
  return(zList)
}


forecast_simulation <- function(object, newxreg, h, nsim, forc_dates, innov = NULL, custom_slope = NULL, init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
{
  switch(object$spec$model$type,
         "1" = forecast_aaa_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, init_states = init_states, innov_type = innov_type, sigma_scale = sigma_scale, ...),
         "2" = forecast_mmm_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, init_states = init_states, innov_type = innov_type, sigma_scale = sigma_scale, ...),
         "3" = forecast_mam_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, init_states = init_states, innov_type = innov_type, sigma_scale = sigma_scale, ...),
         "4" = forecast_powermam_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, init_states = init_states, innov_type = innov_type, sigma_scale = sigma_scale, ...))
}

# forecast_simulation_redraw <- function(object, fcast_dist, newxreg, h, nsim, drop_na, drop_negative, forc_dates, innov = NULL, custom_slope = NULL, init_states = NULL, innov_type = "q", sigma_scale = NULL, ...)
# {
#   nsim_act <- nrow(fcast_dist)
#   resim_iter <- 0
#   max_resim_iter <- 10
#
#   while (nsim_act < nsim && resim_iter < max_resim_iter) {
#     resim_iter <- resim_iter + 1
#
#     resim_ratio <- nsim / max(nsim_act,1)
#     n_resim <- ceiling(resim_ratio * (nsim - nsim_act + 100))
#
#     new_fcast_dist <- forecast_simulation(object = object, newxreg = newxreg, h = h, nsim = n_resim, forc_dates = forc_dates, innov = innov, custom_slope = custom_slope, init_states = init_states, innov_type = innov_type, sigma_scale = sigma_scale, ...)$distribution
#     new_fcast_dist <- forecast_sanitycheck(new_fcast_dist, h, nsim, drop_na, drop_negative, verbose = FALSE)
#     fcast_dist <- rbind(fcast_dist,new_fcast_dist)
#     nsim_act <- nrow(fcast_dist)
#   }
#   if (resim_iter >= max_resim_iter) {
#     warning(paste0("\nmaximum number of resampling iterations reached. Forecast distribution might have less than ", nsim, " draws."))
#   }
#   if (nsim_act > nsim) {
#     # trim excess
#     fcast_dist <- fcast_dist[1:nsim,]
#   }
#   return(fcast_dist)
# }
#
#
# forecast_sanitycheck <- function(fcast_dist, h, nsim, drop_na, drop_negative, verbose = TRUE)
# {
#   if (drop_na) {
#     good_inds <- (rowSums(is.na(fcast_dist) | is.nan(fcast_dist) | is.infinite(fcast_dist)) == 0)
#     if (any(!good_inds)) {
#       if (verbose) {
#         warning("\nNA/NaN/Inf values detected (and removed) from the simulated forecast distribution")
#       }
#       fcast_dist <- fcast_dist[good_inds,]
#       if (length(fcast_dist) < 10*h) {
#         if (verbose) {
#           warning("\nToo many bad draws detected in the simulated forecast distribution; setting forecast matrix to zero")
#         }
#         fcast_dist <- matrix(0, nsim, h)
#       }
#     }
#   }
#   if (drop_negative) {
#     good_inds <- (rowSums( fcast_dist < 0 ) == 0)
#     fcast_dist <- fcast_dist[good_inds,]
#     if (length(fcast_dist) < 10*h) {
#       if (verbose) {
#         warning("\nToo many negative draws detected in the simulated forecast distribution; setting forecast matrix to zero")
#       }
#       fcast_dist <- matrix(0, nsim, h)
#     }
#   }
#   return(fcast_dist)
# }

# mean calculations

forecast_backtransform <- function(fcast_dist, transform)
{
  if (!is.null(transform)) {
    f1 <- transform$inverse(as.numeric(fcast_dist), transform$lambda)
    f1 <- matrix(f1, ncol = ncol(fcast_dist), nrow = nrow(fcast_dist), byrow = FALSE)
    colnames(f1) <- colnames(fcast_dist)
    mean.forecast <- colMeans(f1)
  } else {
    mean.forecast <- colMeans(fcast_dist)
    f1 <- fcast_dist
  }
  return(list(mean = mean.forecast, dist = f1))
}

# class2 and class3 taken from Hyndman's forecast package
# based on sections 6.3 of "Forecasting with Exponential Smoothing" book

analytic_moments.class1 <- function(mod, h = 1, newxreg  = NULL, init_states = NULL)
{
  m <- mod$spec$seasonal$frequency
  if (is.null(init_states)) {
    last_state <- tail(mod$model$states, 1)
  } else {
    last_state <- init_states
  }
  last_state <- as.numeric(last_state)
  par <- coef(mod)
  if (mod$spec$xreg$include_xreg & !is.null(newxreg)) {
    X <- mod$spec$xreg$xreg
    cf <- mod$model$setup$parmatrix
    rnames <- rownames(cf)
    xbeta <- cf[grepl("^rho",rnames),1]
    xreg <- as.numeric(coredata(newxreg) %*% xbeta)
  } else {
    xreg <- rep(0, h)
  }
  sigma2 <- sd(residuals(mod, raw = T))^2
  p <- length(last_state)
  H <- matrix(c(1, rep(0, p - 1)), nrow = 1)
  if (substr(mod$spec$model$model,3,3) == "A") {
    H[1, p] <- 1
  }
  if (substr(mod$spec$model$model,2,2) == "A") {
    if (mod$spec$model$include_damped) {
      H[1, 2] <- par["phi"]
    } else {
      H[1, 2] <- 1
    }
  }
  F <- matrix(0, p, p)
  F[1, 1] <- 1
  if (substr(mod$spec$model$model,2,2) == "A") {
    if (mod$spec$model$include_damped) {
      F[1, 2] <- F[2, 2] <- par["phi"]
    } else {
      F[1, 2] <- F[2, 2] <- 1
    }
  }
  if (substr(mod$spec$model$model,3,3) == "A") {
    F[p - m + 1, p] <- 1
    F[(p - m + 2):p, (p - m + 1):(p - 1)] <- diag(m - 1)
  }
  G <- matrix(0, nrow = p, ncol = 1)
  G[1, 1] <- par["alpha"]
  if (substr(mod$spec$model$model,2,2) == "A") {
    G[2, 1] <- par["beta"]
  }
  if (substr(mod$spec$model$model,3,3) == "A") {
    G[3, 1] <- par["gamma"]
  }
  mu <- numeric(h)
  Fj <- diag(p)
  cj <- numeric(h - 1)
  if (h > 1) {
    for (i in 1:(h - 1))
    {
      mu[i] <- H %*% Fj %*% last_state + xreg[i]
      cj[i] <- H %*% Fj %*% G
      Fj <- Fj %*% F
    }
    cj2 <- cumsum(cj ^ 2)
    var <- sigma2 * c(1, 1 + cj2)
  }
  else {
    var <- sigma2
  }
  mu[h] <- H %*% Fj %*% last_state

  return(list(mu = mu, var = var, cj = cj))

}



analytic_moments.class2 <- function(mod, h = 1, newxreg  = NULL, init_states = NULL)
{
  tmp <- analytic_moments.class2(mod, h = h, newxreg = newxreg, init_states = init_states)
  theta <- numeric(h)
  sigma <- sd(residuals(mod, raw = T))

  theta[1] <- tmp$mu[1] ^ 2
  if (h > 1) {
    for (j in 2:h)
      theta[j] <- tmp$mu[j] ^ 2 + sigma^2 * sum(tmp$cj[1:(j - 1)] ^ 2 * theta[(j - 1):1])
  }
  var <- (1 + sigma^2) * theta - tmp$mu ^ 2
  return(list(mu = tmp$mu, var = var))
}
# need to rewrite this using matrix notation for inclusion of the regressors
# and the array notation

# also need to think about the centered seasonality issue

analytic_moments.class3 <- function(mod, h = 1, newxreg  = NULL, init_states = NULL)
{
  m <- mod$spec$seasonal$frequency
  if (is.null(init_states)) {
    last_state <- tail(mod$model$states, 1)
  } else {
    last_state <- init_states
  }
  p <- length(last_state)
  slope <- ifelse(substr(mod$spec$model$model,2,2) != "N", 1, 0)
  H1 <- matrix(rep(1, 1 + slope), nrow = 1)
  H2 <- matrix(c(rep(0, m - 1), 1), nrow = 1)
  par <- coef(mod)
  sigma2 <- sd(residuals(mod, raw = T))^2
  if (mod$spec$xreg$include_xreg & !is.null(newxreg)) {
    return(analytic_moments.class3x(mod, h = h, newxreg  = newxreg, init_states = init_states))
  }
  if (mod$spec$model$include_damped) {
    phi <- par["phi"]
  } else {
    phi <- 1
  }
  if (slope == 0) {
    F1 <- 1
    G1 <- par["alpha"]
  } else {
    F1 <- rbind(c(1, 1), c(0, 1))
    G1 <- rbind(c(par["alpha"], par["alpha"]), c(par["beta"], par["beta"]))
  }
  F2 <- rbind(c(rep(0, m - 1), 1), cbind(diag(m - 1), rep(0, m - 1)))

  G2 <- matrix(0, m, m)
  G2[1, m] <- par["gamma"]
  Mh <- matrix(last_state[1:(p - m)]) %*% matrix(last_state[(p - m + 1):p], nrow = 1)
  Vh <- matrix(0, length(Mh), length(Mh))
  H21 <- H2 %x% H1
  F21 <- F2 %x% F1
  G21 <- G2 %x% G1
  K <- (G2 %x% F1) + (F2 %x% G1)
  mu <- var <- numeric(h)
  for (i in 1:h)
  {
    mu[i] <- H1 %*% Mh %*% t(H2)
    mu[i] <- mu[i]
    var[i] <- (1 + sigma2) * H21 %*% Vh %*% t(H21) + sigma2 * mu[i] ^ 2
    vecMh <- c(Mh)
    Vh <- F21 %*% Vh %*% t(F21) + sigma2 * (F21 %*% Vh %*% t(G21) + G21 %*% Vh %*% t(F21) + K %*% (Vh + vecMh %*% t(vecMh)) %*% t(K) + sigma2 * G21 %*% (3 * Vh + 2 * vecMh %*% t(vecMh)) %*% t(G21))
    Mh <- F1 %*% Mh %*% t(F2) + G1 %*% Mh %*% t(G2) * sigma2
  }
  return(list(mu = mu, var = var))
}

analytic_moments.class3x <- function(mod, h = 1, newxreg  = NULL, init_states = NULL)
{
  m <- mod$spec$seasonal$frequency
  if (is.null(init_states)) {
    last_state <- tail(mod$model$states, 1)
  } else {
    last_state <- init_states
  }
  last_state <- as.numeric(last_state)
  p <- length(last_state)
  slope <- ifelse(substr(mod$spec$model$model,2,2) != "N", 1, 0)
  par <- coef(mod)
  sigma2 <- sd(residuals(mod, raw = T))^2
  if (mod$spec$xreg$include_xreg & !is.null(newxreg)) {
    cf <- mod$model$setup$parmatrix
    rnames <- rownames(cf)
    xbeta <- cf[grepl("^rho",rnames),1]
    xreg <- as.numeric(coredata(newxreg) %*% xbeta)
  } else {
    xreg <- rep(0, h)
  }
  if (mod$spec$model$include_damped) {
    phi <- cumsum(par["phi"]^(1:h))
  } else {
    phi <- 1:h
  }
  if (mod$spec$model$include_trend) {
    beta <- par["beta"]
    b_n <- last_state[2]
  } else {
    beta <- 0
    b_n <- 0
  }
  s_n <- rev(last_state[(p - m + 1):p])
  gamma <- par["gamma"]
  hplus <- ( ((1:h) - 1) %% m ) + 1
  hm <- ((1:h) - 1)/m
  l_n <- last_state[1]
  alpha <- par["alpha"]
  cj <-  theta <- mu <- var <- rep(0, h)
  for (i in 1:h) {
    nonseasonal_mean <- l_n + phi[i] * b_n + xreg[i]
    mu[i] <- nonseasonal_mean * s_n[hplus[i]]
    cj[i] <- alpha + beta * phi[i]
    if (i == 1) {
      theta[i] <- nonseasonal_mean^2
    } else {
      theta[i] <- nonseasonal_mean^2 + sigma2 * sum((cj[1:(i - 1)]^2) * (theta[i - (1:(i - 1))]))
    }
    var[i] <- s_n[hplus[i]]^2 * (theta[i] * (1 + sigma2) * (1 + gamma^2 * sigma2)^hm[i] - nonseasonal_mean^2)
  }
  return(list(mu = mu, var = var))
}

# ToDO: incomplete -> need to add xreg
tsmoments.tsets.estimate <- function(object, h = 1,  newxreg = NULL, init_states = NULL, ...)
{
  type <- object$model$setup$type
  # no analytical solution or approximation for power MAM
  if (type == 4) return(NULL)
  type <- object$spec$model$class
  # no analytical solution for class 4+ models (MMN, MMM)
  if (type > 3) {
    return(NULL)
  }
  # some checks
  if (!is.null(init_states)) {
    n_required <- ncol(object$model$states)
    if (length(as.numeric(init_states)) != n_required) stop("\nlength of init_states not equal to number of states in model.")
  }
  if (!is.null(newxreg)) {
    include_xreg <- object$spec$model$include_xreg
    if (!include_xreg) {
      warning("\nnewxreg provided when there are no regressors in the model.")
      newxreg <- NULL
    } else {
      n_required <- NCOL(object$spec$xreg$xreg)
      if (NCOL(newxreg) != n_required) stop("\nnumber of columns of newxreg not equal to xreg in model/")
      if (NROW(newxreg) != h) stop("\nnumber of rows of newxreg not equal to horizon (h).")
    }
  }
  out <- switch(as.character(type),
                "1"  = analytic_moments.class1(object, h = h, newxreg = newxreg, init_states = init_states),
                "2" = analytic_moments.class2(object, h = h, newxreg = newxreg, init_states = init_states),
                "3" = analytic_moments.class3(object, h = h, newxreg = newxreg, init_states = init_states))
  return(out)
}
