predict.tsets.estimate <- function(object, h = 12, newxreg = NULL, nsim = 1000, drop_na = TRUE, drop_negative = FALSE, redraw = TRUE, forc_dates = NULL, 
                                   innov = NULL, custom_slope = NULL, ...)
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
  if (!is.null(innov)) {
    requiredn <- h * nsim
    if (length(innov) != requiredn) {
      stop("\nlength of innov must be nsim x h")
    }
    # check that the innovations are uniform samples (from a copula)
    if (any(innov < 0 | innov > 1 )) {
      stop("\ninnov must be >0 and <1 (uniform samples)")
    }
    if (any(innov == 0)) innov[which(innov == 0)] <- 1e-12
    if (any(innov == 1)) innov[which(innov == 1)] <- (1 - 1e-12)
    innov <- matrix(innov, h, nsim)
  }
  # run simulation
  zList <- forecast_simulation(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov = innov, custom_slope = custom_slope, ...)
  # check for bad values in the forecast distribution and resimulate (if necessary)
  if (is.null(innov)) {
    fcast_dist <- forecast_sanitycheck(zList$distribution, h, nsim, drop_na, drop_negative)
    if (redraw && (nrow(fcast_dist) < nsim)) {
      fcast_dist <- forecast_simulation_redraw(object, fcast_dist, newxreg, h, nsim, drop_na, drop_negative, forc_dates = forc_dates, innov = innov,  ...)
    }    
  } else {
    fcast_dist <- zList$distribution
  }
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


forecast_aaa_cpp <- function(object, newxreg = NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, ...)
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
  stateinit <- tail(object$model$states, 1)
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
    E <- t(qnorm(innov, mean = 0, sd = coefficient["sigma"]))
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
  Level <- out$Level[,-1,drop = FALSE]
  Slope <- out$Slope[,-1,drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,-1,drop = FALSE]
  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)
  return(zList)
}

forecast_mmm_cpp <- function(object, newxreg = NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, ...)
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
  stateinit <- tail(object$model$states,1)
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
    E <- t(matrix(qtruncnorm(innov, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim))
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
  #
  #
  Y <- out$Simulated[,-1,drop = FALSE]
  Level <- out$Level[,-1,drop = FALSE]
  Slope <- out$Slope[,-1,drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,-1,drop = FALSE]
  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)
  return(zList)
}


forecast_mam_cpp <- function(object, newxreg=NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, ...)
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
  stateinit <- tail(object$model$states,1)
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
    E <- t(matrix(qtruncnorm(innov, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim))
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
  Level <- out$Level[,-1,drop = FALSE]
  Slope <- out$Slope[,-1,drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,-1,drop = FALSE]
  xregf <- xregf[-1]
  E <- E[,-1,drop = FALSE]
  zList <- wrap_forecast_output(object, Y, Level, Slope, Seasonal, E, xregf, model, forc_dates)

  return(zList)
}

forecast_powermam_cpp <- function(object, newxreg = NULL, h = 12, nsim = 1000, forc_dates = NULL, innov = NULL, custom_slope = NULL, ...)
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
  stateinit <- tail(object$model$states,1)
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
    E <- t(matrix(qtruncnorm(innov, a = -1, mean = 0, sd = coefficient["sigma"]), h, nsim))
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
  #
  Y <- out$Simulated[,-1,drop = FALSE]
  Level <- out$Level[,-1,drop = FALSE]
  Slope <- out$Slope[,-1,drop = FALSE]
  Seasonal <- out$Seasonal
  Seasonal <- t(Seasonal[,frequency,])
  Seasonal <- Seasonal[,-1,drop = FALSE]
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
  date_class <- attr(object$spec$target$sampling, "date_class")
  dist_classes <- c("tsets.distribution", "tsmodel.distribution")
  class(Level) <- dist_classes
  attr(Level, "date_class") <- date_class
  class(Y)     <- dist_classes
  attr(Y, "date_class") <- date_class
  class(E)     <- dist_classes
  attr(E, "date_class") <- date_class
  if (model[1] == 0) {
    Slope <- NULL
  } else {
    class(Slope) <- dist_classes
    attr(Slope, "date_class") <- date_class
  }
  if (model[2] == 0) {
    Seasonal <- NULL
  } else{
    class(Seasonal) <- dist_classes
    attr(Seasonal, "date_class") <- date_class
  }
  if (object$model$setup$include_xreg == 0) {
    xregf <- NULL
  }
  decomposition <- list(Level = Level, Slope = Slope, Seasonal = Seasonal, X = xregf, Error = E, dates = as.character(forc_dates))
  zList <- list(distribution = Y, original_series = zoo(object$spec$target$y_orig, object$spec$target$index), h = model[3], 
                spec = object$spec, decomposition = decomposition)
  class(zList) <- dist_classes
  return(zList)
}


forecast_simulation <- function(object, newxreg, h, nsim, forc_dates, innov = NULL, custom_slope = NULL, ...)
{
  switch(object$spec$model$type,
         "1" = forecast_aaa_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, ...),
         "2" = forecast_mmm_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, ...),
         "3" = forecast_mam_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, ...),
         "4" = forecast_powermam_cpp(object = object, newxreg = newxreg, h = h, nsim = nsim, forc_dates = forc_dates, innov, custom_slope = custom_slope, ...))
}

forecast_simulation_redraw <- function(object, fcast_dist, newxreg, h, nsim, drop_na, drop_negative, forc_dates, innov = NULL, custom_slope = NULL, ...)
{
  nsim_act <- nrow(fcast_dist)
  resim_iter <- 0
  max_resim_iter <- 10

  while (nsim_act < nsim && resim_iter < max_resim_iter) {
    resim_iter <- resim_iter + 1

    resim_ratio <- nsim / max(nsim_act,1)
    n_resim <- ceiling(resim_ratio * (nsim - nsim_act + 100))

    new_fcast_dist <- forecast_simulation(object = object, newxreg = newxreg, h = h, nsim = n_resim, forc_dates = forc_dates, innov = innov, custom_slope = custom_slope, ...)$distribution
    new_fcast_dist <- forecast_sanitycheck(new_fcast_dist, h, nsim, drop_na, drop_negative, verbose = FALSE)
    fcast_dist <- rbind(fcast_dist,new_fcast_dist)
    nsim_act <- nrow(fcast_dist)
  }
  if (resim_iter >= max_resim_iter) {
    warning(paste0("\nmaximum number of resampling iterations reached. Forecast distribution might have less than ", nsim, " draws."))
  }
  if (nsim_act > nsim) {
    # trim excess
    fcast_dist <- fcast_dist[1:nsim,]
  }
  return(fcast_dist)
}


forecast_sanitycheck <- function(fcast_dist, h, nsim, drop_na, drop_negative, verbose = TRUE)
{
  if (drop_na) {
    good_inds <- (rowSums(is.na(fcast_dist) | is.nan(fcast_dist) | is.infinite(fcast_dist)) == 0)
    if (any(!good_inds)) {
      if (verbose) {
        warning("\nNA/NaN/Inf values detected (and removed) from the simulated forecast distribution")
      }
      fcast_dist <- fcast_dist[good_inds,]
      if (length(fcast_dist) < 10*h) {
        if (verbose) {
          warning("\nToo many bad draws detected in the simulated forecast distribution; setting forecast matrix to zero")
        }
        fcast_dist <- matrix(0, nsim, h)
      }
    }
  }
  if (drop_negative) {
    good_inds <- (rowSums( fcast_dist < 0 ) == 0)
    fcast_dist <- fcast_dist[good_inds,]
    if (length(fcast_dist) < 10*h) {
      if (verbose) {
        warning("\nToo many negative draws detected in the simulated forecast distribution; setting forecast matrix to zero")
      }
      fcast_dist <- matrix(0, nsim, h)
    }
  }
  return(fcast_dist)
}

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

analytic_moments.aaa <- function(mod, h = 1)
{
  m <- mod$spec$seasonal$frequency
  l_n <- tail(mod$model$states, 1)[,"Level"]
  alpha <- coef(mod)["alpha"]
  if (mod$spec$model$include_trend == 1) {
    b_n <- tail(mod$model$states, 1)[,"Trend"]
    beta <- coef(mod)["beta"]
  } else {
    b_n <- 0
    beta <- 0
  }
  if (mod$spec$model$include_damped == 1) {
    phi <- coef(mod)["phi"]
  } else {
    phi <- 1
  }
  if (mod$spec$model$include_seasonal == 1) {
    s_n <- rev(tail(mod$model$states, 1)[,grepl("S[0-9]",colnames(mod$model$states))])
    gamma <- coef(mod)["gamma"]
  } else {
    s_n <- NULL
    gamma <- 0
  }
  sigma <- sd(residuals(mod, raw = T))
  hplus <- ( ((1:h) - 1) %% m ) + 1
  mu <- rep(0, h)
  v <- rep(0, h)
  c <- rep(0, h)
  for (i in 1:h) {
    if (i %% m == 0) d <- 1 else d <- 0
    c[i] <- alpha + beta*sum(phi^(1:i)) + gamma * d
    mu[i] <- l_n + sum(phi^(1:i)) * b_n
    if (!is.null(s_n)) {
      mu[i] <- mu[i] + s_n[hplus[i]]
    }
    if (i == 1) {
      v[i] <- sigma^2
    } else {
      v[i] <- sigma^2 * (1 + sum(c[1:(i - 1)]^2))
    }
  }
  return(list(mu = mu, var = v, cj = c))
}

analytic_moments.mam <- function(mod, h = 1)
{
  m <- mod$spec$seasonal$frequency
  l_n <- tail(mod$model$states, 1)[,"Level"]
  alpha <- coef(mod)["alpha"]
  if (mod$spec$model$include_trend == 1) {
    b_n <- tail(mod$model$states, 1)[,"Trend"]
    beta <- coef(mod)["beta"]
  } else {
    b_n <- 0
    beta <- 0
  }
  if (mod$spec$model$include_damped == 1) {
    phi <- coef(mod)["phi"]
  } else {
    phi <- 1
  }
  if (mod$spec$model$include_seasonal == 1) {
    s_n <- rev(tail(mod$model$states, 1)[,grepl("S[0-9]",colnames(mod$model$states))])
    gamma <- coef(mod)["gamma"]
  } else {
    s_n <- NULL
    gamma <- 0
  }
  sigma <- sd(residuals(mod, raw = T))
  hplus <- ( ((1:h) - 1) %% m ) + 1
  mu <- rep(0, h)
  v <- rep(0, h)
  c <- rep(0, h)
  for (i in 1:h) {
    if (i %% m == 0) d <- 1 else d <- 0
    c[i] <- alpha + beta*sum(phi^(1:i)) + gamma * d
    mu[i] <- l_n + sum(phi^(1:i)) * b_n
    if (!is.null(s_n)) {
      mu[i] <- mu[i] + s_n[hplus[i]]
    }
    if (i == 1) {
      v[i] <- sigma^2
    } else {
      v[i] <- (1+sigma^2) * (1 + sum(c[1:(i - 1)]^2))
    }
  }
  return(list(mu = mu, var = v, cj = c))
}

analytic_moments.mmm <- function(mod, h = 1)
{
  m <- mod$spec$seasonal$frequency
  l_n <- tail(mod$model$states, 1)[,"Level"]
  alpha <- coef(mod)["alpha"]
  if (mod$spec$model$include_trend == 1) {
    b_n <- tail(mod$model$states, 1)[,"Trend"]
    beta <- coef(mod)["beta"]
  } else {
    b_n <- 0
    beta <- 0
  }
  if (mod$spec$model$include_damped == 1) {
    phi <- coef(mod)["phi"]
  } else {
    phi <- 1
  }
  if (mod$spec$model$include_seasonal == 1) {
    s_n <- rev(tail(mod$model$states, 1)[,grepl("S[0-9]",colnames(mod$model$states))])
    gamma <- coef(mod)["gamma"]
  } else {
    s_n <- NULL
    gamma <- 0
  }
  sigma <- sd(residuals(mod, raw = T))
  hplus <- ( ((1:h) - 1) %% m ) + 1
  mu <- rep(0, h)
  v <- rep(0, h)
  c <- rep(0, h)
  for (i in 1:h) {
    if (i %% m == 0) d <- 1 else d <- 0
    c[i] <- alpha + beta*sum(phi^(1:i)) + gamma * d
    mu[i] <- l_n + sum(phi^(1:i)) * b_n
    if (!is.null(s_n)) {
      mu[i] <- mu[i] + s_n[hplus[i]]
    }
    if (i == 1) {
      v[i] <- sigma^2
    } else {
      v[i] <- sigma^2 * (1 + sum(c[1:(i - 1)]^2))
    }
  }
  return(list(mu = mu, var = v, cj = c))
}

