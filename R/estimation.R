#' Model Estimation
#'
#' @description Estimates a model given a specification object using
#' maximum likelihood.
#' @param object an object of class \dQuote{tsets.spec}.
#' @param solver one of either \dQuote{solnp}, \dQuote{nlminb} or \dQuote{optim}.
#' The latter uses the L-BFGS-B algorithm from the lbfgsb3c package. For option
#' \dQuote{autodiff}, valid solvers are \dQuote{nlminb} and \dQuote{nloptr}.
#' @param control solver control parameters.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsetsad package.
#' @param ... only additional argument which can be passed when choosing autodiff
#' is that of \dQuote{use_hessian} which tells the solver to make use of second
#' derivatives.
#' @details The maximum likelihood estimation uses bound constraints for some of
#' the parameters as described in the online book of tsmodels. Additionally, for
#' parameters which are constrained to be less than another parameter
#' (e.g. beta<alpha), a simple barrier approach is adopted which adjusts the
#' previous likelihood value upwards by some fixed percentage of that value
#' during the minimization. The observation variance is not directly estimated
#' but instead concentrated out of the likelihood. When autodiff is TRUE with
#' the nloptr solver, the constraints and their jacobian are explicitly used,
#' whilst a soft barrier constraint is used in the case of the nlminb solver.
#' @return An object of class \dQuote{tsets.estimate}
#' @aliases estimate
#' @method estimate tsets.spec
#' @rdname estimate
#' @export
#'
#'
estimate.tsets.spec <- function(object, solver = "nlminb", control = list(trace = 0), autodiff = TRUE, ...)
{
  # create an environment for the soft barrier solver
  ets_env <- new.env(hash = TRUE)
  assign("ets_llh", 1, envir = ets_env)
  setup <- object$model
  pars <- setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"init"]
  setup$parnames <- names(pars)
  setup$data <- as.numeric(object$target$y)
  setup$good <- rep(1, length(setup$data))
  setup$good[which(is.na(setup$data))] <- 0
  setup$good <- c(1,setup$good)
  setup$frequency <- object$seasonal$frequency
  setup$xreg <- object$xreg$xreg
  setup$k <- ncol(object$xreg$xreg)
  setup$ets_env <- ets_env
  lfun <- tsets_ll_fun(setup$type)
  ffun <- tsets_filter_fun(setup$type)
  setup$estimation <- 1
  # find parameter scaling
  if (any(grepl("rho",names(pars)))) {
    ix <- which(grepl("rho",names(pars)))
    xr <- setup$xreg
    yr <- setup$data
    pscale <- sapply(1:length(ix), function(i) max(yr, na.rm = TRUE)/max(xr[,i]))
    if (any(!is.finite(pscale))) {
      pscale[which(!is.finite(pscale))] <- 1
    }
    parscale <- rep(1,length(pars))
    parscale[ix] <- pscale
  } else {
    parscale <- rep(1, length(pars))
  }
  tic <- Sys.time()
    # run solver
  setup$debug <- FALSE
  hess <- NULL
  if (autodiff) {
    if (!solver %in% c("nlminb","nloptr")) {
      solver <- "nlminb"
      warning("For auto_diff, only nlminb and nloptr allowed. Setting solver to nlminb.\n")
    }
    opt_res <- estimate_ad(object, solver = solver, control = control, ...)
    pars <- opt_res$pars
    lik <- opt_res$llh
    hess <- opt_res$hess
    opt_res <- opt_res$solver_out
  } else {
    if (solver == "nlminb") {
      opt_res <- nlminb(start = pars, objective = lfun, lower = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"lower"],
                        upper = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"upper"], setup = setup,
                        scale = 1/parscale, control = control)
      pars <- opt_res$par
      lik <- opt_res$objective
      #
      if (opt_res$convergence != 0) {
        # try again
        run_count <- 0
        cont_ <- TRUE
        while (cont_) {
          run_count <- run_count + 1
          opt_res <- nlminb(start = pars, objective = lfun, lower = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"lower"],
                            upper = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"upper"], setup = setup,
                            scale = 1/parscale, control = control)
          pars <- opt_res$par
          lik <- opt_res$objective
          if (opt_res$convergence == 0 || run_count > 19) {
            cont_ <- FALSE
          }
        }
        pars <- opt_res$par
        lik <- opt_res$objective
      }
      if (opt_res$convergence != 0) warnings("\nnlminb indicates non-successful convergence...")
    } else if (solver == "solnp") {
      opt_res <- solnp(pars = pars, fun = lfun, LB = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"lower"],
                       UB = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"upper"], setup = setup, control = control)
      pars <- opt_res$pars
      lik <- opt_res$values[length(opt_res$values)]
      if (opt_res$convergence != 0) warnings("\nsolnp indicates non successful convergence...")
    } else if (solver == "optim") {
      control$parscale <- parscale
      mult_trend <- (setup$type == 2)
      mult_seasonality <- (setup$type == 2) || (setup$type == 3) || (setup$type == 4)
      setup$impose_bounds <- TRUE
      inp_pars_trans <- pars_estim_pre_trans(pars,mult_trend,mult_seasonality, setup$parmatrix)
      opt_res <- optim(par = inp_pars_trans, fn = lfun, setup = setup, method = "Nelder-Mead", control = control)
      #
      setup$impose_bounds <- FALSE
      pars_inp <- pars_estim_inv_trans(opt_res$par,mult_trend,mult_seasonality, setup$parmatrix)
      opt_res <- optim(par = pars_inp, fn = lfun, lower = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"lower"],
                       upper = setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"upper"], setup = setup,
                       method = "L-BFGS-B", control = control, hessian = F)
      pars <- opt_res$par
      lik <- opt_res$value
      if (opt_res$convergence != 0) warnings("\noptim indicates non-successful convergence...")
    } else {
      stop("\nsolver must be either solnp, nlminb or optim")
    }
  }
  if (object$target$scaler != 1 & object$model$type == 1) {
    pnames <- names(pars)
    pars["l0"] <- pars["l0"] * object$target$scaler
    if (any(pnames == "b0")) pars["b0"] <- pars["b0"] * object$target$scaler
    if (any(grepl("rho",pnames))) pars[grepl("rho",pnames)] <- pars[grepl("rho",pnames)] * object$target$scaler
    if (any(grepl("[s][0-9]",pnames))) pars[grepl("[s][0-9]",pnames)] <- pars[grepl("[s][0-9]",pnames)] * object$target$scaler
    setup$data <- setup$data * object$target$scaler
    object$target$scaler <- 1
    lik <- lfun(pars, setup)
  }
  setup$estimation <- 0
  setup$parmatrix[which(setup$parmatrix[,"estimate"] == 1),"init"] <- pars
  if (setup$include_xreg == 1) {
    x <- c(0, setup$xreg %*% setup$parmatrix[paste0("rho",1:ncol(setup$xreg)),1])
  } else {
    if (setup$type == 2) {
      x <- rep(0, length(setup$data) + 1)
    } else {
      x <- rep(0, length(setup$data) + 1)
    }
  }

  if (setup$type == 4) {
    filt <- ffun(setup$data, alpha = setup$parmatrix["alpha","init"], beta = setup$parmatrix["beta","init"],
                 gamma = setup$parmatrix["gamma","init"], phi = setup$parmatrix["phi","init"],
                 l0 = setup$parmatrix["l0","init"], b0 = setup$parmatrix["b0","init"],
                 theta = setup$parmatrix["theta","init"], delta = setup$parmatrix["delta","init"],
                 s0 = setup$parmatrix[paste0("s",0:(setup$frequency - 2)),"init"],
                 frequency = setup$frequency, x = x, setup = setup)
  } else {
    filt <- ffun(setup$data, alpha = setup$parmatrix["alpha","init"], beta = setup$parmatrix["beta","init"],
                 gamma = setup$parmatrix["gamma","init"], phi = setup$parmatrix["phi","init"],
                 l0 = setup$parmatrix["l0","init"], b0 = setup$parmatrix["b0","init"],
                 s0 = setup$parmatrix[paste0("s",0:(setup$frequency - 2)),"init"],
                 frequency = setup$frequency, x = x, setup = setup)
  }

  setup$parmatrix["sigma","init"] <- sd(filt$residuals[which(setup$good[-1] == 1)])
  filt$setup <- setup
  filt$loglik <- lik
  opt_res$timing <- difftime(Sys.time(), tic, units = "mins")
  obj <- list(model = filt, spec = object, opt = opt_res, hess = hess)
  class(obj) <- c("tsets.estimate","tsmodel.estimate")
  return(obj)
}

tsets_ll_fun <- function(type)
{
  fun <- switch(as.character(type),
               "1" = tsets_ll_aaa,
               "2" = tsets_ll_mmm,
               "3" = tsets_ll_mam,
               "4" = tsets_ll_powermam)
  return(fun)
}

tsets_filter_fun <- function(type)
{
  fun <- switch(as.character(type),
               "1" = tsets_filter_aaa,
               "2" = tsets_filter_mmm,
               "3" = tsets_filter_mam,
               "4" = tsets_filter_powermam)
  return(fun)
}


# Log-likelihood functions

tsets_ll_aaa <- function(pars, setup)
{
  ets_env <- setup$ets_env
  names(pars) <- setup$parnames
  if (!is.null(setup$impose_bounds)) {
    if (setup$impose_bounds) {
      pars <- pars_estim_inv_trans(pars,FALSE,FALSE,setup$parmatrix)
    }
  }
  parmatrix <- setup$parmatrix
  parmatrix[which(parmatrix[,"estimate"] == 1),"init"] <- pars
  l0 <- parmatrix["l0","init"]
  b0 <- parmatrix["b0","init"]
  s0 <- parmatrix[paste0("s",0:(setup$frequency - 2)),"init"]
  alpha <- parmatrix["alpha","init"]
  beta <- parmatrix["beta","init"]
  gamma <- parmatrix["gamma","init"]
  phi <- parmatrix["phi","init"]
  rho <- parmatrix[paste0("rho",1:setup$k),"init"]

  xreg <- setup$xreg
  if (setup$include_xreg == 1) {
    x <- as.numeric(setup$xreg %*% rho)
    x <- c(0, x)
  } else {
    x <- rep(0, length(setup$data) + 1)
  }

  # soft barrier
  if (setup$estimation == 1) {
    if (beta  >= (alpha - 1e-6)) {
      return(get("ets_llh", ets_env) + 0.2 * (abs(get("ets_llh", ets_env))))
    }
    if (gamma >= (1 - alpha - 1e-6)) {
      return(get("ets_llh", ets_env) + 0.2 * (abs(get("ets_llh", ets_env))))
    }
  }
  pars <- as.numeric(parmatrix[c("l0","b0","alpha","beta","gamma","phi"),1])
  model <- c(setup$include_trend, setup$include_seasonal, setup$frequency, NROW(setup$data) + 1, setup$normalized_seasonality)
  y <- c(0, as.numeric(setup$data))

  if (!setup$debug) {
    f <- filter_aaa(model_ = model, y_ = y, pars_ = pars, s0_ = s0, x_ = x, good_ = setup$good)
    if (setup$estimation == 1) {
      if (any(is.na(f$Error))) {
        return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
      }
    }
    n <- sum(setup$good) - 1
    Error <- f$Error[which(setup$good == 1)]
    out <- n * log(sum(Error[-1]^2))
    if (setup$estimation == 1) assign("ets_llh", out, envir = ets_env)
  } else {
    f <- tsets_filter_aaa(setup$data, alpha = alpha, beta = beta, gamma = gamma, phi = phi, l0 = l0, b0 = b0, s0 = s0, frequency = setup$frequency, x = x, setup = setup)
    if (setup$estimation == 1) {
      if (any(is.na(f$residuals))) {
        return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
      }
    }
    n <- length(setup$data[which(setup$good == 1)])
    Error <- f$residuals[which(setup$good == 1)]
    out <- n * log(sum(Error^2))
    if (setup$estimation == 1) {
      assign("ets_llh", out, envir = ets_env)
    }
  }
  return(out)
}

tsets_ll_mmm <- function(pars, setup)
{
  ets_env <- setup$ets_env
  names(pars) <- setup$parnames

  if (!is.null(setup$impose_bounds)) {
    if (setup$impose_bounds) {
      pars <- pars_estim_inv_trans(pars,TRUE,TRUE,setup$parmatrix)
    }
  }
  parmatrix <- setup$parmatrix
  parmatrix[which(parmatrix[,"estimate"] == 1),"init"] <- pars
  alpha <- parmatrix["alpha","init"]
  l0 <- parmatrix["l0","init"]
  theta <- parmatrix["theta","init"]
  b0 <- parmatrix["b0","init"]
  beta <- parmatrix["beta","init"]
  phi <- parmatrix["phi","init"]
  s0 <- parmatrix[paste0("s",0:(setup$frequency - 2)),"init"]
  gamma <- parmatrix["gamma","init"]
  rho <- parmatrix[paste0("rho",1:setup$k),"init"]

  xreg <- setup$xreg
  if (setup$include_xreg == 1) {
    x <- as.numeric(setup$xreg %*% rho)
    x <- c(0, x)
  } else {
    x <- rep(0, length(setup$data) + 1)
  }
  pars <- as.numeric(parmatrix[c("l0","b0","alpha","beta","gamma","phi"),1])
  model <- c(setup$include_trend, setup$include_seasonal, setup$frequency, NROW(setup$data) + 1, setup$normalized_seasonality)
  y <- c(0, as.numeric(setup$data))

  if (!setup$debug) {
    f <- filter_mmm(model_ = model, y_ = y, pars_ = pars, s0_ = s0, x_ = x, good_ = setup$good)
    if (setup$estimation == 1) {
      if (any(is.na(f$Error))) {
        return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
      }
    }
    n <- sum(setup$good) - 1
    Error <- f$Error[which(setup$good == 1)]
    Filt <- f$Filtered[which(setup$good == 1)]
    out <- n * log(sum(Error[-1]^2)) + 2 * sum(log(abs(Filt[-1])))
    if (setup$estimation == 1) {
      assign("ets_llh", out, envir = ets_env)
    }
  } else {
    f <- tsets_filter_mmm(setup$data, alpha, beta, gamma, phi, l0, b0, s0, frequency = setup$frequency, x = x, setup = setup)
    if (setup$estimation == 1) {
      if (any(is.na(f$residuals))) {
        return(get("ets_llh", ets_env) + 0.3 * (abs(get("ets_llh", ets_env))))
      }
    }
    n <- length(setup$data)
    out <- n * log(sum(f$residuals^2)) + 2 * sum(log(abs(f$filtered)))

    if (setup$estimation == 1) {
      assign("ets_llh", out, envir = ets_env)
    }
  }
  return(out)
}

tsets_ll_mam <- function(pars, setup)
{
  ets_env <- setup$ets_env
  names(pars) <- setup$parnames

  if (!is.null(setup$impose_bounds)) {
    if (setup$impose_bounds) {
      pars <- pars_estim_inv_trans(pars,FALSE,TRUE,setup$parmatrix)
    }
  }
  parmatrix <- setup$parmatrix
  parmatrix[which(parmatrix[,"estimate"] == 1),"init"] <- pars
  alpha <- parmatrix["alpha","init"]
  l0 <- parmatrix["l0","init"]
  theta <- parmatrix["theta","init"]
  b0 <- parmatrix["b0","init"]
  beta <- parmatrix["beta","init"]
  phi <- parmatrix["phi","init"]
  s0 <- parmatrix[paste0("s",0:(setup$frequency - 2)),"init"]
  gamma <- parmatrix["gamma","init"]
  rho <- parmatrix[paste0("rho",1:setup$k),"init"]
  sig <- parmatrix["sigma",1]

  xreg <- setup$xreg
  if (setup$include_xreg == 1) {
    x <- as.numeric(setup$xreg %*% rho)
    x <- c(0, x)
  } else {
    x <- rep(0, length(setup$data) + 1)
  }

  if (setup$estimation == 1) {
    if (beta >= (alpha - 1e-6)) {
      return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
    }
  }

  pars <- as.numeric(parmatrix[c("l0","b0","alpha","beta","gamma","phi"),1])
  model <- c(setup$include_trend, setup$include_seasonal, setup$frequency, NROW(setup$data) + 1, setup$normalized_seasonality)
  y <- c(0, as.numeric(setup$data))

  if (!setup$debug) {
    f <- filter_mam(model_ = model, y_ = y, pars_ = pars, s0_ = s0, x_ = x, good_ = setup$good)
    if (setup$estimation == 1) {
      if (any(is.na(f$Error))) {
        return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
      }
    }
    n <- sum(setup$good) - 1
    Error <- f$Error[which(setup$good == 1)]
    Filt <- f$Filtered[which(setup$good == 1)]
    out <- n * log(sum(Error[-1]^2)) + 2 * sum(log(abs(Filt[-1])))
    if (setup$estimation == 1) {
      assign("ets_llh", out, envir = ets_env)
    }
  } else {
    f <- tsets_filter_mam(setup$data, alpha, beta, gamma, phi, l0, b0, s0, frequency = setup$frequency, x = x, setup = setup)

    if (setup$estimation == 1) {
      if (any(is.na(f$residuals))) {
        return(get("ets_llh", ets_env) + 0.3 * (abs(get("ets_llh", ets_env))))
      }
    }
    n <- length(setup$data)
    out <- n * log(sum(f$residuals^2)) + 2 * sum(log(abs(f$filtered)))

    if (setup$estimation == 1) {
      assign("ets_llh", out, envir = ets_env)
    }
  }
  return(out)
}

tsets_ll_powermam <- function(pars, setup)
{
  ets_env <- setup$ets_env

  if (!is.null(setup$impose_bounds)) {
    if (setup$impose_bounds) {
      pars <- pars_estim_inv_trans(pars,FALSE,TRUE,setup$parmatrix)
    }
  }

  parmatrix <- setup$parmatrix
  parmatrix[which(parmatrix[,"estimate"] == 1),"init"] <- pars

  alpha <- parmatrix["alpha","init"]
  beta <- parmatrix["beta","init"]
  rho <- parmatrix[paste0("rho",1:setup$k),"init"]
  s0 <- parmatrix[paste0("s",0:(setup$frequency - 2)),"init"]

  xreg <- setup$xreg
  if (setup$include_xreg == 1) {
    x <- as.numeric(setup$xreg %*% rho)
    x <- c(0, x)
  } else {
    x <- rep(0, length(setup$data) + 1)
  }
  if (setup$estimation == 1) {
    if (beta >= (alpha - 1e-6)) return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
  }

  pars <- as.numeric(parmatrix[c("l0","b0","alpha","beta","gamma","phi","theta","delta"),1])
  model <- c(setup$include_trend, setup$include_seasonal, setup$frequency, NROW(setup$data) + 1, setup$normalized_seasonality)
  y <- c(0, as.numeric(setup$data))

  f <- filter_powermam(model_ = model, y_ = y, pars_ = pars, s0_ = s0, x_ = x, good_ = setup$good)
  if (setup$estimation == 1) {
    if (any(is.na(f$Error))) return(get("ets_llh", ets_env) + 0.1 * (abs(get("ets_llh", ets_env))))
  }
  n <- sum(setup$good) - 1
  Error <- f$Error[which(setup$good == 1)]
  Filt <- f$Fpower[which(setup$good == 1)]
  out <- n * log(sum(Error[-1]^2)) + 2 * sum(log(abs(Filt[-1])))
  if (setup$estimation == 1) assign("ets_llh", out, envir = ets_env)
  return(out)
}


# Functions to impose parameter bounds for optim Nelder-Mead step (contribution from Keith O'Hara)

gen_logit_trans <- function(val_inp,lower_bound,upper_bound)
{
  eps_dbl <- 0

  if (val_inp == lower_bound) {
    par_val_out <- -Inf
  } else if (val_inp == upper_bound) {
    par_val_out <- Inf
  } else {
    par_val_out <- log(val_inp - lower_bound + eps_dbl) - log(upper_bound - val_inp + eps_dbl)
  }

  return(par_val_out)
}

gen_inv_logit_trans <- function(val_trans_inp,lower_bound,upper_bound)
{
  eps_dbl <- 0

  if (val_trans_inp == -Inf) {
    par_val_out <- lower_bound
  } else if (val_trans_inp == Inf) {
    par_val_out <- upper_bound
  } else {
    par_val_out <- (lower_bound + eps_dbl + (upper_bound - eps_dbl)*exp(val_trans_inp)) / (1.0 + exp(val_trans_inp))
  }

  if (is.nan(par_val_out)) {
    if (val_trans_inp < 0) {
      par_val_out <- lower_bound
    } else {
      par_val_out <- upper_bound
    }
  }

  return(par_val_out)
}

pars_estim_pre_trans <- function(pars, mult_trend, mult_season, parmatrix)
{
  pars_names <- names(pars)

  if ("alpha" %in% pars_names) {
    orig_alph <- pars["alpha"]
    pars["alpha"] <- gen_logit_trans(pars["alpha"], parmatrix["alpha","lower"], parmatrix["alpha","upper"])
  }

  if ("beta" %in% pars_names) {
    if (mult_trend) {
      pars["beta"] <- gen_logit_trans(pars["beta"], parmatrix["beta","lower"], parmatrix["beta","upper"])
    } else {
      apar <- parmatrix["alpha","init"]
      pars["beta"] <- gen_logit_trans(pars["beta"], 0, apar - 1e-6)
    }
  }

  if ("gamma" %in% pars_names) {
    if (mult_season) {
      pars["gamma"] <- gen_logit_trans(pars["gamma"], parmatrix["gamma","lower"], parmatrix["gamma","upper"])
    } else {
      apar <- parmatrix["alpha","init"]
      pars["gamma"] <- gen_logit_trans(pars["gamma"], 0, 1 - apar - 1e-6)
    }
  }

  if ("phi" %in% pars_names) {
    pars["phi"] <- gen_logit_trans(pars["phi"], parmatrix["phi","lower"], parmatrix["phi","upper"])
  }

  return(pars)
}

pars_estim_inv_trans <- function(pars, mult_trend, mult_season, parmatrix)
{
  pars_names <- names(pars)

  if ("alpha" %in% pars_names) {
    pars["alpha"] <- gen_inv_logit_trans(pars["alpha"], parmatrix["alpha","lower"], parmatrix["alpha","upper"])
  }
  if ("beta" %in% pars_names) {
    if (mult_trend) {
      pars["beta"] <- gen_inv_logit_trans(pars["beta"], parmatrix["beta","lower"], parmatrix["beta","upper"])
    } else {
      apar <- parmatrix["alpha","init"]
      pars["beta"] <- gen_inv_logit_trans(pars["beta"], 0, apar)
    }
  }

  if ("gamma" %in% pars_names) {
    if (mult_season) {
      pars["gamma"] <- gen_inv_logit_trans(pars["gamma"], parmatrix["gamma","lower"], parmatrix["gamma","upper"])
    } else {
      apar <- parmatrix["alpha","init"]
      pars["gamma"] <- gen_inv_logit_trans(pars["gamma"], 0, 1 - apar - 1e-6)
    }
  }

  if ("phi" %in% pars_names) {
    pars["phi"] <- gen_inv_logit_trans(pars["phi"], parmatrix["phi","lower"], parmatrix["phi","upper"])
  }

  return(pars)
}
#################################################################################
# From Rob Hyndman's forecast package (not yet used).
admissible_region <- function(alpha, beta = NULL, gamma = NULL, phi = 1, m = 1) {
  if (phi < 0 || phi > 1 + 1e-8) {
    return(FALSE)
  }
  if (is.null(gamma)) {
    if (alpha < 1 - 1 / phi || alpha > 1 + 1 / phi) {
      return(FALSE)
    }
    if (!is.null(beta)) {
      if (beta < alpha * (phi - 1) || beta > (1 + phi) * (2 - alpha)) {
        return(FALSE)
      }
    }
  } else if (m > 1) {
    if (is.null(beta)) {
      beta <- 0
    }
    if (gamma < max(1 - 1 / phi - alpha, 0) || gamma > 1 + 1 / phi - alpha) {
      return(FALSE)
    }
    if (alpha < 1 - 1 / phi - gamma * (1 - m + phi + phi * m) / (2 * phi * m)) {
      return(FALSE)
    }
    if (beta < -(1 - phi) * (gamma / m + alpha)) {
      return(FALSE)
    }
    P <- c(phi * (1 - alpha - gamma), alpha + beta - alpha * phi + gamma - 1, rep(alpha + beta - alpha * phi, m - 2), (alpha + beta - phi), 1)
    roots <- polyroot(P)
    if (max(abs(roots)) > 1 + 1e-10) {
      return(FALSE)
    }
  }
  # pass
  return(TRUE)
}
#################################################################################
