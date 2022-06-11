#' Model Simulation
#'
#' @description Simulates paths from an ETS model with optional user specified values.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param nsim the number of paths per complete set of time steps (h).
#' @param seed an object specifying if and how the random number generator should be
#' initialized (\sQuote{seeded}). Either NULL or an integer that will be used in
#' a call to set.seed before simulating the response vectors. If set, the value
#' is saved as the "seed" attribute of the returned value. The default, NULL
#' will not change the random generator state, and return .Random.seed as the
#' \dQuote{seed} attribute in the returned object.
#' @param h the number of time steps to simulate paths for. If this is NULL,
#' it will use the same number of periods as in the original series.
#' @param newxreg an optional matrix of regressors to use for the simulation
#' if xreg was used in the estimation. If NULL and the estimated object had
#' regressors, and h was also set to NULL, then the original regressors will
#' be used.
#' @param sim_dates an optional vector of simulation dates equal to h. If NULL
#' will use the implied periodicity of the data to generate a regular sequence
#' of dates after the first available date in the data.
#' @param bootstrap whether to bootstrap the innovations from the estimated
#' object by re-sampling from the empirical distribution.
#' @param innov an optional vector of uniform innovations which will be
#' translated to regular innovations using the appropriate distribution quantile
#' function and model standard deviation. The length of this vector should be
#' equal to nsim x horizon.
#' @param sigma_scale an optional scalar which will scale the standard deviation
#' of the innovations (useful for profiling under different assumptions). This
#' applies to all approaches (parametric, bootstrap and custom innovations supplied
#' values).
#' @param pars an optional named vector of model coefficients which override the
#' estimated coefficients. No checking is currently performed on the adequacy of
#' these coefficients.
#' @param ... not currently used.
#' @details It is required that an initial object of class \dQuote{tsets.estimate}
#' be passed to this function, but the parameters and initial states can be overriden
#' by passing them as named values in the \sQuote{pars} argument. The default is
#' to initialize the states from the seed states used in the estimated object,
#' with h equal to the length of the original series (default for NULL h).
#' Innovations for the simulation can either be parametric (Normal for additive
#' or truncated Normal for multiplicative error models), based on the estimated
#' residuals (bootstrap argument) or use supplied set of uniform random numbers
#' (which are then translated into Normal or truncated Normal using standard
#' deviation equal to the model sigma and optionally scaled by sigma_scale).
#' @return An object of class \dQuote{tsets.simulate} with slots for the simulated
#' series and states.
#' @aliases simulate
#' @method simulate tsets.estimate
#' @rdname simulate
#' @export
#'
#'
simulate.tsets.estimate  = function(object, nsim = 1, seed = NULL, h = NULL, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, innov = NULL,
                                    sigma_scale = 1, pars = coef(object), ...)
{
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  if (is.null(h)) h <- length(object$spec$target$y)
  frequency <- ifelse(object$model$setup$include_seasonal == 1, object$spec$target$frequency, 4)
  ipars <- object$model$setup$parmatrix
  if (!is.null(pars)) {
    if (length(pars[na.omit(match(rownames(ipars),names(pars)))]) == 0) stop("\nno pars names matched with model parameters")
    ipars[names(pars),"init"] <- pars[na.omit(match(rownames(ipars),names(pars)))]
    object$model$setup$parmatrix <- ipars
  }

  if (object$spec$xreg$include_xreg == 1) {
    n_xreg <- ncol(object$spec$xreg$xreg)
    if (ncol(newxreg) != n_xreg) stop("\nnewxreg does not have the same number of regressors as used in the model.")
    h <- nrow(newxreg)
    if (is.xts(newxreg)) {
      sim_dates <- index(newxreg)
    } else {
      if (is.null(sim_dates)) {
        sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
      } else {
        if (length(sim_dates) != h) stop("\nsim_dates not equal to nrow xreg")
      }
    }
    if (!all.equal(sort(colnames(newxreg)),sort(colnames(object$spec$xreg$xreg)))) {
      stop("\nnewxreg colnames not the same as xreg used in model")
    }
    newxreg <- newxreg[,colnames(object$spec$xreg$xreg)]
    newxreg <- coredata(newxreg)
    rhovector <- ipars[grepl("[rho][0-9]",rownames(ipars)),1,drop = FALSE]
    x <- newxreg %*% rhovector
  } else {
    if (is.null(sim_dates)) {
      sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
    } else {
      if (length(sim_dates) != h) stop("\nsim_dates not equal to nrow xreg")
    }
    x <- matrix(rep(0, h), ncol = 1)
  }
  date_class <- attr(object$spec$target$sampling, "date_class")

  x <- rbind(matrix(0, ncol = ncol(x), nrow = 1), x)
  custom_flag <- 0
  model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim, custom_flag)

  sigma_res <- object$model$setup$parmatrix["sigma",1]
  type <- object$spec$model$type

  if (bootstrap) {
    res <- as.numeric(na.omit(object$model$residuals))
    res <- res * sigma_scale
    E <- matrix(sample(res, (h * nsim), replace = TRUE), ncol = h, nrow = nsim)
  } else {
    if (type == 1) {
      if (is.null(innov)) {
        E <- matrix(rnorm(h * nsim, 0, sigma_res * sigma_scale), ncol = h, nrow = nsim)
      } else {
        if (length(innov) != (h * nsim)) {
          stop("\nlength innov must be nsim x h")
        }
        if (any(innov) == 0) {
          innov[which(innov == 0)] <- 1e-12
        }
        if ( any(innov == 1)) {
          innov[which(innov == 1)] <- 1 - 1e-12
        }
        innov <- matrix(innov, h, nsim)
        E <- qnorm(innov, sd = sigma_res * sigma_scale)
      }
    } else {
      if (is.null(innov)) {
        E <- matrix(pmax(-1, rnorm(h * nsim, 0, sigma_res * sigma_scale)), ncol = h, nrow = nsim)
      } else {
        if (length(innov) != (h * nsim)) {
          stop("\nlength innov must be nsim x h")
        }
        if (any(innov) == 0) {
          innov[which(innov == 0)] <- 1e-12
        }
        if ( any(innov == 1)) {
          innov[which(innov == 1)] <- 1 - 1e-12
        }
        innov <- matrix(innov, h, nsim)
        E <- qtruncnorm(innov, mean = 0, sd = sigma_res * sigma_scale, a = -1)
      }
    }
  }
  E <- cbind(matrix(0, nrow = nsim, ncol = 1), E)

  # extract parameters
  alpha <- ipars["alpha",1]
  beta <- ipars["beta",1]
  gamma <- ipars["gamma",1]
  phi <- ipars["phi",1]
  delta <- ipars["delta",1]
  theta <- ipars["theta",1]
  svector <- ipars[grepl("s[0-9]",rownames(ipars)),1]
  if (object$spec$model$error == "Additive") {
    svector <- c(svector, -sum(svector))
  } else {
    svector <- c(svector, frequency - sum(svector))
  }
  l0 <- ipars["l0",1]
  b0 <- ipars["b0",1]
  if (type == 4) {
    pars <- rep(0, 8)
    pars[1] <- l0
    pars[2] <- b0
    pars[3] <- alpha
    pars[4] <- beta
    pars[5] <- gamma
    pars[6] <- phi
    pars[7] <- theta
    pars[8] <- delta
  } else {
    pars <- rep(0, 6)
    pars[1] <- l0
    pars[2] <- b0
    pars[3] <- alpha
    pars[4] <- beta
    pars[5] <- gamma
    pars[6] <- phi
  }
  out <- switch(type,
               "1" = simulate_aaa(model_ = model, e_ = E, pars_ = pars, s0_ = svector, x_ = x, slope_overide_ = matrix(1)),
               "2" = simulate_mmm(model_ = model, e_ = E, pars_ = pars, s0_ = svector, x_ = x, slope_overide_ = matrix(1)),
               "3" = simulate_mam(model_ = model, e_ = E, pars_ = pars, s0_ = svector, x_ = x, slope_overide_ = matrix(1)),
               "4" = simulate_powermam(model_ = model, e_ = E, pars_ = pars, s0_ = svector, x_ = x, slope_overide_ = matrix(1)))
  Y <- out$Simulated[,-1]
  Level <- out$Level[,-1]
  Slope <- out$Slope[,-1]
  Seasonal <- out$Seasonal[-1,,]
  Seasonal <- t(Seasonal[, frequency, ])
  x <- x[-1]
  E <- E[,-1]
  sclasses <- c("tsets.distribution","tsmodel.distribution")
  if (type == 1) {
      f1 <- forecast_backtransform(Y, object$spec$transform)
      Y <- f1$dist
  }
  colnames(Y) <- as.character(sim_dates)
  colnames(Level) <- as.character(sim_dates)
  colnames(Slope) <- as.character(sim_dates)
  colnames(Seasonal) <- as.character(sim_dates)
  names(x) <- as.character(sim_dates)
  class(Level) <- sclasses
  attr(Level, "date_class") <- date_class
  class(Y) <- sclasses
  attr(Y, "date_class") <- date_class
  class(E) <- sclasses
  attr(E, "date_class") <- date_class
  if (model[1] == 0) {
    Slope <- NULL
  } else{
    class(Slope) <- sclasses
    attr(Slope, "date_class") <- date_class
  }
  if (model[2] == 0) {
    Seasonal <- NULL
  } else {
    class(Seasonal) <- sclasses
    attr(Seasonal, "date_class") <- date_class
  }
  if (object$spec$model$include_xreg == 0) {
    x <- NULL
  }
  zList <- list(Simulated = Y, Level = Level, Slope = Slope, Seasonal = Seasonal, X = x, Error = E, dates = as.character(sim_dates),
               seed = RNGstate, pars = ipars[ipars[,"required"] == 1, "init"], sigma = ipars["sigma","init"], sigma_scale = sigma_scale)
  class(zList) <- c("tsets.simulate")
  return(zList)
}



ets_sample <- function(model = "AAA", power = FALSE, damped = FALSE, h = 100, frequency = 12, start_date = as.Date("1990-01-01"),
                       alpha = 0.06, beta = 0.01, phi = 1, gamma = 0.01, delta = 1, theta = 1,
                       rho = NULL, xreg = NULL, normalized_seasonality = TRUE, sigma = 1,
                       transformation = NULL, lambda = NULL, lower = 0, upper = 1,
                       seed_states = c("l0" = 1, "b0" = 0.05), innov = NULL, seed = NULL, sampling = "months")
{
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  if (is.null(h)) {
    h <- 100
  } else {
    h <- max(1, h)
  }
  if (is.null(start_date)) {
    start_date <- as.Date("1990-01-01")
    date_class <- "Date"
  } else {
    date_class <- is(start_date)
    if (!any(date_class %in% c("Date","POSIXct"))) stop("start_date must be either Date or POSIXct")
  }
  s_dates <- future_dates(start = start_date, sampling, n = h)
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

  # check init_states

  if (substr(model, 1, 1) == "A" & power) {
    warning("\npower model not available with additive error...setting to FALSE")
    power <- FALSE
  }
  if (substr(model,2,2) == "M" & power) {
    warning("\npower model not available with multiplicative trend...setting to FALSE")
    power <- FALSE
  }
  error_type <- ifelse(substr(model,1,1) == "A", "Additive","Multiplicative")
  model <- validate_model(model)
  type <- model_type(model, power)
  damped <- as.logical(damped)
  include_trend <- ifelse(substr(model,2,2) != "N", 1, 0)
  include_seasonal <- ifelse(substr(model,3,3) != "N", 1, 0)
  include_damped <- ifelse(include_trend == 1 && damped, 1, 0)
  if (!is.null(frequency)) {
    frequency <- max(1, as.integer(frequency))
  } else {
    frequency <- 4
  }
  if (include_seasonal == 1) {
    if (is.null(frequency)) stop("\nseasonal model requested but frequency is NULL.")
    if (frequency == 1) stop("\nseasonal model requested but frequency is 1.")
  } else {
    frequency <- 4
  }
  s_vector_names <-  paste0("s",0:(frequency - 2))

  if (!is.null(xreg)) {
    xreg <- coredata(xreg)
    if (nrow(xreg) != h) stop("\nxreg rows must equal h")
    if (is.null(rho)) stop("\nxreg present but rho (coefficients of xreg) is NULL")
    if (length(rho) != ncol(xreg)) stop("\nlength of rho not equal to cols of xreg")
    include_xreg <- 1
    x <- xreg %*% rho
  } else {
    include_xreg <- 0
    rho <- NULL
    x <- matrix(0, ncol = 1, nrow = h)
  }
  x <- rbind(matrix(0, ncol = ncol(x), nrow = 1), x)

  # check stability of coefficients

  if (is.null(alpha[1])) alpha <- 0.06
  if (is.null(beta[1])) beta <- 0.01
  if (is.null(phi[1])) phi <- 1
  if (is.null(gamma[1])) gamma <- 0.01
  if (is.null(delta[1])) delta <- 1
  if (is.null(theta[1])) theta <- 1
  if (is.null(seed_states["l0"])) {
    l0 <- 10
  } else {
    l0 <- seed_states["l0"]
  }
  if (is.null(seed_states["b0"])) {
    b0 <- 0.5
  } else {
    b0 <- seed_states["b0"]

  }
  if (is.null(sigma[1])) sigma <- l0 * 0.1
  if (include_seasonal) {
    if (any(grepl(pattern = "s[0-9]", names(seed_states)))) {
      s_vector <- unlist(seed_states[grepl(pattern = "s[0-9]", names(seed_states))])
      user_s_vector_names <- sort(names(s_vector))
      if (!all(user_s_vector_names %in% s_vector_names) | length(user_s_vector_names) != length(s_vector_names)) {
        cat("\nseed_states for seasonal do not match required names :", s_vector_names)
        stop("\nExiting")
      } else {
        s_vector <- s_vector[s_vector_names]
      }
    } else {
      if (error_type == "Additive") {
        s_vector <- fourier_series(as.Date(1:(10 * length(s_vector_names)), origin = "1970-01-01"), period = frequency, K = length(s_vector_names)/2)
        s_vector <- (s_vector %*% rnorm(ncol(s_vector), mean = 0.005 * l0, sd = 0.05 * l0))[1:length(s_vector_names),1]
        names(s_vector) <- s_vector_names
      } else {

      }
    }
  } else {
    if (error_type == "Additive") {
      s_vector <- rep(0, length(s_vector_names))
      names(s_vector) <- s_vector_names
    } else {
      s_vector <- rep(1, length(s_vector_names))
      names(s_vector) <- s_vector_names
    }
  }

  p_matrix <- data.table("parameters" = c("l0","b0",s_vector_names,"alpha","beta","phi","gamma","delta","theta","sigma"),
                         values = c(seed_states["l0"],seed_states["b0"], s_vector, alpha[1], beta[1], phi[1], gamma[1], delta[1], theta[1], sigma[1]))

  check_table <- check_parameters_simulation(p_matrix, model)
  if (any(!check_table$`>lb`) | any(!check_table$`<ub`) | any(!na.omit(check_table$condition_pass))) {
    warning("\nparameters violate stability conditions: ")
    print(check_table)
  }
  custom_flag <- 0
  model <- c(include_trend, include_seasonal, h, frequency, normalized_seasonality, 1, custom_flag)
  if (type == 1) {
    if (is.null(innov)) {
      E <- matrix(rnorm(h * 1, 0, sigma), ncol = h, nrow = 1)
    } else {
      if (length(innov) != (h * 1)) {
        stop("\nlength innov must be h")
      }
      if (any(innov) == 0) {
        innov[which(innov == 0)] <- 1e-12
      }
      if ( any(innov == 1)) {
        innov[which(innov == 1)] <- 1 - 1e-12
      }
      innov <- matrix(innov, h, 1)
      E <- qnorm(innov, sd = sigma)
    }
  } else {
    if (is.null(innov)) {
      E <- matrix(pmax(-1, rnorm(h * 1, 0, sigma)), ncol = h, nrow = 1)
    } else {
      if (length(innov) != (h * 1)) {
        stop("\nlength innov must be h")
      }
      if (any(innov) == 0) {
        innov[which(innov == 0)] <- 1e-12
      }
      if ( any(innov == 1)) {
        innov[which(innov == 1)] <- 1 - 1e-12
      }
      innov <- matrix(innov, h, 1)
      E <- qtruncnorm(innov, mean = 0, sd = sigma, a = -1)
    }
  }
  E <- cbind(matrix(0, nrow = 1, ncol = 1),  E)

  # extract parameters
  parameters <- NULL
  alpha <- p_matrix[parameters == "alpha"]$values
  beta <- p_matrix[parameters == "beta"]$values
  gamma <- p_matrix[parameters == "gamma"]$values
  phi <- p_matrix[parameters == "phi"]$values
  delta <- p_matrix[parameters == "delta"]$values
  theta <- p_matrix[parameters == "theta"]$values
  if (error_type == "Additive") {
    s_vector <- c(s_vector, -sum(s_vector))
  } else {
    s_vector <- c(s_vector, frequency - sum(s_vector))
  }
  l0 <- p_matrix[parameters == "l0"]$values
  b0 <- p_matrix[parameters == "b0"]$values
  if (type == 4) {
    pars <- rep(0, 8)
    pars[1] <- l0
    pars[2] <- b0
    pars[3] <- alpha
    pars[4] <- beta
    pars[5] <- gamma
    pars[6] <- phi
    pars[7] <- theta
    pars[8] <- delta
  } else {
    pars <- rep(0, 6)
    pars[1] <- l0
    pars[2] <- b0
    pars[3] <- alpha
    pars[4] <- beta
    pars[5] <- gamma
    pars[6] <- phi
  }
  out <- switch(type,
                "1" = simulate_aaa(model_ = model, e_ = E, pars_ = pars, s0_ = s_vector, x_ = x, slope_overide_ = matrix(1)),
                "2" = simulate_mmm(model_ = model, e_ = E, pars_ = pars, s0_ = s_vector, x_ = x, slope_overide_ = matrix(1)),
                "3" = simulate_mam(model_ = model, e_ = E, pars_ = pars, s0_ = s_vector, x_ = x, slope_overide_ = matrix(1)),
                "4" = simulate_powermam(model_ = model, e_ = E, pars_ = pars, s0_ = s_vector, x_ = x, slope_overide_ = matrix(1)))
  Y <- out$Simulated[,-1]
  Level <- out$Level[,-1]
  Slope <- out$Slope[,-1]
  Seasonal <- out$Seasonal[-1,,,drop = FALSE]
  Seasonal <- Seasonal[, frequency, ]
  x <- x[-1]
  E <- E[,-1]
  s_names <- c("Simulated","Level")
  if (type == 1) {
    if (!is.null(transformation)) {
        if (transformation[1] == "box-cox" & is.na(lambda)) stop("\nlambda cannot be NA in simulation")
        transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
        Y <- transform$inverse(Y, lambda = lambda)
    }
  }
  if (model[1] == 0) {
    Slope <- NULL
  } else {
    s_names <- c(s_names, "Slope")
  }
  if (model[2] == 0) {
    Seasonal <- NULL
  } else {
    s_names <- c(s_names, "Seasonal")
  }
  if (include_xreg == 0) {
    x <- NULL
  } else {
    s_names <- c(s_names, "X")
  }
  s_names <- c(s_names,"Error")
  sim <- cbind(Y, Level, Slope, Seasonal, x, E)
  colnames(sim) <- s_names
  sim <- xts(sim, s_dates)
  return(sim)
}



check_parameters_simulation <- function(pars, model) {
  error_type <- ifelse(substr(model,1,1) == "A", "Additive","Multiplicative")
  cf_names <- pars$parameters
  cf <- pars$values
  names(cf) <- cf_names
  if (error_type == "Additive") {
    if (any(cf_names == "beta")) {
      condition_slope <- cf["beta"]  <= (cf["alpha"] - 0.01)
    } else {
      condition_slope <- NA
    }
    if (any(cf_names == "gamma")) {
      condition_seasonal <- cf["gamma"]  <= (1 - cf["alpha"] - 0.01)
    } else {
      condition_seasonal <- NA
    }
  }
  lb_check <- cf[c("alpha","beta","gamma","phi","theta","delta")] >= 1e-12
  ub_check <- cf[c("alpha","beta","gamma","phi","theta","delta")] <= 1

  if (error_type == "Additive") {
    condition_table <- data.table(coef = c("alpha","beta","gamma","phi","theta","delta"), value = round(as.numeric(cf[c("alpha","beta","gamma","phi","theta","delta")]),4),
                                  ">lb" = lb_check, "<ub" = ub_check, "condition" = c("NA"," < alpha"," < (1 - alpha)","NA","NA","NA"), "condition_pass" = c(NA, condition_slope, condition_seasonal,NA, NA, NA))
  } else {
    condition_table <- data.table(coef = c("alpha","beta","gamma","phi","theta","delta"), value = c(round(as.numeric(cf[c("alpha","beta","gamma","phi","theta","delta")]),4)),
                                  ">lb" = c(lb_check), "<ub" = c(ub_check), "condition_pass" = TRUE)
  }
  return(condition_table)
}
