validate_model <- function(model = "AAA")
{
  valid_models <- c("AAA","AAN","ANN","ANA","MMM","MMN","MNM","MNN","MAM","MAN","Theta")
  if (any(!model %in% valid_models)) {
    stop("\nNot a valid model.")
  } else {
    return(model)
  }
}

model_type <- function(model = "AAA", power = FALSE)
{
  if (substr(model,1,1) == "A") {
    type <- 1
  } else if (substr(model,1,1) == "M") {
    type <- 2
  } else {
    type <- 2
  }
  if (substr(model,1,2) == "MA") {
    type <- 3
    if (power) type <- 4
  }
  return(type)
}

init_pars <- function(alpha = NULL, beta = NULL, gamma = NULL, phi = NULL, trend_type = "A", season_type = "M", damped = FALSE, lower = c(rep(1e-4, 3), 0.8), upper = rep(0.99,4), frequency = 12)
{
  if (any(lower > upper)) {
    stop("Inconsistent parameter boundaries")
  }
  # alpha
  if (is.null(alpha)) {
    alpha <- lower[1] + 0.2 * (upper[1] - lower[1]) / frequency
    if (alpha > 1 || alpha < 0) {
      alpha <- lower[1] + 2e-3
    }
    par <- c(alpha = alpha)
  } else {
    par <- numeric(0)
  }
  # beta
  if (trend_type != "N" && is.null(beta)) {
    # Ensure beta < alpha
    upper[2] <- min(upper[2], alpha)
    beta <- lower[2] + 0.1 * (upper[2] - lower[2])
    if (beta < 0 || beta > alpha) {
      beta <- alpha - 1e-3
    }
    par <- c(par, beta = beta)
  }
  # gamma
  if (season_type != "N" && is.null(gamma)) {
    # Ensure gamma < 1-alpha
    upper[3] <- min(upper[3], 1 - alpha)
    gamma <- lower[3] + 0.05 * (upper[3] - lower[3])
    if (gamma < 0 || gamma > 1 - alpha) {
      gamma <- 1 - alpha - 1e-3
    }
    par <- c(par, gamma = gamma)
  }

  # phi
  if (damped && is.null(phi)) {
    phi <- lower[4] + .99 * (upper[4] - lower[4])
    if (phi < 0 || phi > 1) {
      phi <- upper[4] - 1e-3
    }
    par <- c(par, phi = phi)
  }

  return(par)
}

init_states <- function(y, trend_type = "A", season_type = "M", frequency = 12){
  if (season_type != "N") {
    # Do decomposition
    n <- NROW(y)
    if (n < 4) {
      stop("not enough data for a seasonal model!!!")
    } else if (n < (3 * frequency)) {
      fouriery <- fourier_series(y, period = frequency, K = 1)
      fit <- tslinear(y, trend = TRUE,  xreg = fouriery)
      if (season_type == "A") {
        y_d <- list(seasonal = y - fit$coef[1] - fit$coef[2] * (1:n))
      } else {
        # season_type=="M". Biased method, but we only need a starting point
        y_d <- list(seasonal = y / (fit$coef[1] + fit$coef[2] * (1:n)))
      }
    } else {
      # n is large enough to do a decomposition
      y_d <- decompose(y, type = switch(season_type, A = "additive", M = "multiplicative"))
    }
    # Initial seasonal component
    init_seas <- rev(y_d$seasonal[2:frequency])
    names(init_seas) <- paste("s", 0:(frequency - 2), sep = "")
    # Seasonally adjusted data
    if (season_type == "A") {
      y_sa <- y - y_d$seasonal
    } else {
      # We do not want negative seasonal indexes
      init_seas <- pmax(init_seas, 1e-2) 
      if (sum(init_seas) > frequency) {
        init_seas <- init_seas / sum(init_seas + 1e-2)
      }
      y_sa <- y / pmax(y_d$seasonal, 1e-2)
    }
  } else{
    frequency <- 1
    init_seas <- NULL
    y_sa <- y
  }
  maxn <- min(max(10, 2 * frequency), length(y_sa))
  if (trend_type == "N") {
    l0 <- mean(y_sa[1:maxn])
    b0 <- NULL
  } else {
    # Simple linear regression on seasonally adjusted data
    fit <- lsfit(1:maxn, y_sa[1:maxn])
    if (trend_type == "A") {
      l0 <- fit$coef[1]
      b0 <- fit$coef[2]
      # If error type is "M", then we don't want l0+b0=0,
      # so perturb just in case.
      if (abs(l0 + b0) < 1e-8) {
        l0 <- l0 * (1 + 1e-3)
        b0 <- b0 * (1 - 1e-3)
      }
    } else {
      # First fitted value
      l0 <- fit$coef[1] + fit$coef[2]
      if (abs(l0) < 1e-8) {
        l0 <- 1e-7
      }
      # Ratio of first two fitted values
      b0 <- (fit$coef[1] + 2 * fit$coef[2]) / l0
      # First fitted value divided by b0
      l0 <- l0 / b0
      if (abs(b0) > 1e10) {
        # Avoid infinite slopes
        b0 <- sign(b0) * 1e10
      }
      if (l0 < 1e-8 || b0 < 1e-8) {
        # Simple linear approximation didn't work.
        l0 <- max(y_sa[1], 1e-3)
        b0 <- max(y_sa[2] / y_sa[1], 1e-3)
      }
    }
  }
  names(l0) <- "l0"
  if (!is.null(b0)) {
    names(b0) <- "b0"
  }
  return(c(l0, b0, init_seas))
}

init_model = function(y, model = "AAA", damped = FALSE, power = FALSE, xreg = NULL, frequency = 12, positive = FALSE)
{
  pars <- "alpha"
  lower <- 1e-04
  upper <- 1 - 1e-6
  doestimate <- 1
  required <- 1
  init <- 0
  # Trend
  if (substr(model,2,2) != "N") {
    include_trend <- TRUE
    pars <- c(pars,"beta")
    lower <- c(lower,1e-04)
    upper <- c(upper,1 - 1e-6)
    doestimate <- c(doestimate,1)
    required <- c(required,1)

    if (substr(model,2,2) == "M") {
      init <- c(init,0)
    } else {
      init <- c(init,0)
    }
  } else {
    pars <- c(pars,"beta")
    lower <- c(lower,1e-04)
    upper <- c(upper,1 - 1e-6)
    doestimate <- c(doestimate,0)
    include_trend <- FALSE
    required <- c(required,0)
    
    if (substr(model,2,2) == "M") {
      init <- c(init,0)
    } else {
      init <- c(init,0)
    }
  }

  # Seasonality
  if (substr(model,3,3) != "N") {
    include_season <- TRUE
    pars <- c(pars,"gamma")
    lower <- c(lower,0)
    upper <- c(upper,1 - 1e-6)
    doestimate <- c(doestimate,1)
    required <- c(required, 1)
    if (substr(model,3,3) == "M") {
      init <- c(init,0)
    } else {
      init <- c(init,0)
    }
  } else {
    include_season <- FALSE
    pars <- c(pars,"gamma")
    lower <- c(lower,0)
    upper <- c(upper,1 - 1e-6)
    doestimate <- c(doestimate,0)
    required <- c(required,0)

    if (substr(model,3,3) == "M") {
      init <- c(init,0)
    } else {
      init <- c(init,0)
    }
  }

  # dampening
  if (include_trend & damped) {
    pars <- c(pars,"phi")
    lower <- c(lower,0.8)
    upper <- c(upper,1)
    doestimate <- c(doestimate,1)
    required <- c(required, 1)
    if (substr(model,2,2) == "M") {
      init <- c(init,1)
    } else {
      init <- c(init,1)
    }
  } else {
    pars <- c(pars,"phi")
    lower <- c(lower,0.8)
    upper <- c(upper,1)
    doestimate <- c(doestimate,0)
    required <- c(required, 0)
    if (substr(model,2,2) == "M") {
      init <- c(init,1)
    } else {
      init <- c(init,1)
    }
  }

  # Power models
  if (substr(model,1,2) == "MA" & power) {
    if (substr(model,3,3) != "N") {
      pars <- c(pars,"theta","delta")
      lower <- c(lower,0,0)
      upper <- c(upper,1,1)
      doestimate <- c(doestimate,1,1)
      required <- c(required,1,1)
      init <- c(init,1,1)
    } else {
      pars <- c(pars,"theta","delta")
      lower <- c(lower,0,0)
      upper <- c(upper,1,1)
      doestimate <- c(doestimate,1,0)
      required <- c(required,1,0)
      init <- c(init,1,1)
    }
  } else {
    pars <- c(pars,"theta","delta")
    lower <- c(lower,0,0)
    upper <- c(upper,1,1)
    doestimate <- c(doestimate,0,0)
    required <- c(required,0,0)
    init <- c(init,1,1)
  }

  #
  # Initialization

  # l0
  pars <- c(pars,"l0")
  init <- c(init,0)

  if (substr(model,1,1) == "A") {
    if (positive) {
      lower <- c(lower,0)
    } else {
      lower <- c(lower,-Inf)
    }
    upper <- c(upper,Inf)
    doestimate <- c(doestimate,1)
    required <- c(required, 1)
  } else {
    lower <- c(lower,1e-12)
    upper <- c(upper,Inf)
    doestimate <- c(doestimate,1)
    required <- c(required,1)
  }

  # b0
  if (include_trend) {
    pars <- c(pars, "b0")
    if (substr(model,2,2) == "M") {
      init <- c(init,1)
      lower <- c(lower,0)
      upper <- c(upper,2)
      doestimate <- c(doestimate,1)
      required = c(required,1)
    } else {
      if (substr(model,1,1) == "A") {
        init <- c(init,0)
        lower <- c(lower,-Inf)
        upper <- c(upper,Inf)
        doestimate <- c(doestimate,1)
        required <- c(required, 1)
      } else {
        init <- c(init,0)
        lower <- c(lower, 0)
        upper <- c(upper,Inf)
        doestimate <- c(doestimate,1)
        required <- c(required, 1)
      }
    }
  } else {
    pars <- c(pars,"b0")
    lower <- c(lower,-Inf)
    upper <- c(upper, Inf)
    doestimate <- c(doestimate,0)
    required <- c(required, 0)

    if (substr(model,2,2) == "M") {
      init <- c(init,1)
    } else {
      if (substr(model,1,1) == "M" & substr(model,2,2) == "N") {
        init <- c(init,1)
      } else {
        init <- c(init,0)
      }
    }
  }

  # s0
  if (include_season) {
    pars <- c(pars,paste0("s", 0:(frequency - 2)))
    if (substr(model,3,3) == "M") {
      lower <- c(lower, rep(0, frequency - 1))
      upper <- c(upper, rep(Inf, frequency - 1))
    } else {
      lower <- c(lower, rep(-Inf, frequency - 1))
      upper <- c(upper, rep( Inf, frequency - 1))
    }
    if (frequency <= 52) {
      doestimate <- c(doestimate, rep(1,frequency - 1))
    } else {
      doestimate <- c(doestimate, rep(0,frequency - 1))
    }

    required <- c(required, rep(1,frequency - 1))
    if (substr(model,3,3) == "M") {
      init <- c(init, rep(1, frequency - 1))
    } else {
      init <- c(init, rep(0,frequency - 1))
    }
  } else {
    if (is.null(frequency)) frequency <- 4
    pars <- c(pars, paste0("s", 0:(frequency - 2)))
    lower <- c(lower, rep(NA,frequency - 1))
    upper <- c(upper, rep(NA,frequency - 1))
    doestimate <- c(doestimate, rep(0, frequency - 1))
    required <- c(required, rep(0, frequency - 1))
    if (substr(model,1,1) == "M" | substr(model,3,3) == "M") {
      init <- c(init, rep(1, frequency - 1))
    } else {
      init <- c(init, rep(0, frequency - 1))
    }
  }

  # External regressors
  if (!is.null(xreg)) {
    k <- ncol(xreg)
    pars <- c(pars, paste0("rho",1:k))
    if (substr(model,2,2) == "M" | (substr(model,1,1) == "M" & substr(model,2,2) == "N")) {
      init <- c(init, rep(0,k))
      lower <- c(lower,rep(-Inf,k))
      upper <- c(upper,rep( Inf,k))
      doestimate <- c(doestimate,rep(1,k))
      required <- c(required, rep(1,k))
    } else {
      init <- c(init, rep(0,k))
      lower <- c(lower,rep(-Inf,k))
      upper <- c(upper,rep(Inf,k))
      doestimate <- c(doestimate,rep(1,k))
      required <- c(required, rep(1,k))
    }
  } else {
    k <- 1
    pars <- c(pars, paste0("rho",1:k))
    if (substr(model,2,2) == "M") {
      init <- c(init, rep(0,k))
      lower <- c(lower,rep(-Inf,k))
      upper <- c(upper,rep(Inf,k))
      doestimate <- c(doestimate,rep(0,k))
      required <- c(required, rep(0,k))
    } else {
      init <- c(init, rep(0,k))
      lower <- c(lower,rep(-Inf,k))
      upper <- c(upper,rep(Inf,k))
      doestimate <- c(doestimate,rep(0,k))
      required <- c(required, rep(0,k))
    }
  }

  #
  pars <- c(pars,"sigma")
  init <- c(init, 1)
  lower <- c(lower,1e-12)
  upper <- c(upper,Inf)
  doestimate <- c(doestimate,0)
  required <- c(required, 0)
  parmatrix <- matrix(0, ncol = 6, nrow = length(pars))
  rownames(parmatrix) <- pars
  colnames(parmatrix) <- c("init","lower","upper","estimate","required","fixed")

  # initialize to something
  parmatrix[,"init"] <- init
  ipars <- init_pars(alpha = NULL, beta = NULL, gamma = NULL, phi = NULL, trend_type = substr(model,2,2), 
                    season_type = substr(model,3,3), damped = damped, lower = c(0.2,0.1,0.1,0.5), upper = c(1,1,1,1), 
                    frequency = frequency)

  istart <- init_states(ts(as.numeric(y), frequency = frequency), substr(model,2,2), season_type = substr(model,3,3), frequency = frequency)
  parmatrix[names(ipars),"init"] <- ipars
  parmatrix[names(istart),"init"] <- istart
  parmatrix[pars,"lower"] <- lower
  parmatrix[pars,"upper"] <- upper
  parmatrix[pars,"estimate"] <- doestimate
  parmatrix[pars,"required"] <- required
  parmatrix["sigma","init"] <- 1
  parmatrix[which(parmatrix[,"estimate"] == 0 & parmatrix[,"required"] == 1),"fixed"] <- 1
  parmatrix[which(parmatrix[,"estimate"] == 0 & parmatrix[,"required"] == 0),"init"] <- init[which(parmatrix[,"estimate"] == 0 & parmatrix[,"required"] == 0)]

  return(parmatrix)
}
