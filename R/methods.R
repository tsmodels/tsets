#' Model Fitted Values
#'
#' @description Extract the fitted values from an estimated model.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param ... not currently used.
#' @aliases fitted
#' @method fitted tsets.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.tsets.estimate <- function(object, ...)
{
  if (!is.null(object$spec$transform)) {
    f <- object$spec$transform$inverse(object$model$filtered, object$spec$transform$lambda)
  } else {
    f <- object$model$filtered
  }
  out <- xts(f, object$spec$target$index)
  colnames(out) <- "fitted"
  return(out)
}

#' Model Residuals
#'
#' @description Extract the residual values from an estimated model.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param raw raw residuals are the model based values which for the additive
#' model are on the Box Cox scale, whilst for multiplicative models are equal
#' to actual/fitted - 1.
#' @param h the horizon (steps) ahead residuals required. The default represents
#' the standard residuals whilst for h>1 these are the (1:h)-step ahead in-sample
#' predicted residuals for each time point under fixed coefficients.
#' @param seed a seed value which initializes the simulated predictive distribution
#' from which the h-step ahead forecasts are made in order to calculate the residuals.
#' @param trace whether to show the progress bar for the h-step ahead residuals
#' calculation. The user is expected to have set up appropriate handlers for
#' this using the \dQuote{progressr} package.
#' @param ... not currently used.
#' @details For h>1, this is like performing an in-sample backtest starting at
#' time 1 with fixed coefficients. The purpose of having the matrix of h-step
#' ahead residuals is in order to calculate the 1:h covariance matrix as well
#' as the cross 1:h covariance matrix when ensembling series at multiple horizons.
#' @return An xts vector of the model residuals for h = 1, else a data.table
#' with rows representing the first prediction date and columns the h-ahead
#' forecast residuals.
#' @note The function can use parallel functionality (for h>1) as long as the
#' user has set up a \code{\link[future]{plan}} using the future package.
#' @aliases residuals
#' @method residuals tsets.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.tsets.estimate <- function(object, raw = FALSE, h = 1, seed = NULL, trace = FALSE, ...)
{
  if (h > 1) {
    out <- hresiduals.tsets.estimate(object, h = h, seed = seed, trace = trace, raw = raw, ...)
  } else{
    if (raw) {
      if (object$spec$model$error == "Additive") {
        if (!is.null(object$spec$transform)) {
          out <- xts(object$spec$transform$transform(object$spec$target$y_orig, object$spec$transform$lambda) - object$model$filtered, object$spec$target$index)
        } else {
          out <- xts(object$spec$target$y_orig - object$model$filtered, object$spec$target$index)
        }
      } else {
        out <- xts(object$spec$target$y_orig/object$model$filtered - 1, object$spec$target$index)
      }
    } else {
      out <- xts(object$spec$target$y_orig - as.numeric(fitted(object)), object$spec$target$index)
    }
    colnames(out) <- "residuals"
  }
  return(out)
}

#' Model Log-Likelihood
#'
#' @description Extract the log-likelihood from an estimated model.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param ... not currently used.
#' @return Returns an object of class logLik. This is a number with at least one
#' attribute, "df" (degrees of freedom), giving the number of (estimated)
#' parameters in the model.
#' @aliases logLik
#' @method logLik tsets.estimate
#' @rdname logLik
#' @export
#'
#'
logLik.tsets.estimate <- function(object, ...)
{
  np <- NROW(object$model$setup$parmatrix[which(object$model$setup$parmatrix[,"estimate"] == 1),])
  structure(-0.5 * object$model$loglik, df = np + 1, class = "logLik")
}


#' Performance Metrics
#'
#' @description Performance metrics from an estimated or predicted tsets model.
#' @param object an object of class \dQuote{tsets.estimate} or \dQuote{tsets.predict}
#' @param actual the actual data matched to the dates of the forecasts.
#' @param alpha the coverage level for distributional forecast metrics.
#' @param ... not currently used.
#' @aliases tsmetrics
#' @method tsmetrics tsets.predict
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.tsets.predict = function(object, actual, alpha = 0.1, ...)
{
  n <- NCOL(object$distribution)
  if (NROW(actual) != n) stop("\nactual length not equal to forecast length")
  m_mape <- mape(actual, colMeans(object$distribution))
  m_bias <- bias(actual, colMeans(object$distribution))
  m_mslre <- mslre(actual, colMeans(object$distribution))
  m_mase <- mase(actual, colMeans(object$distribution), object$original_series, frequency = object$spec$target$frequency)
  if (!is.null(alpha)) {
    m_mis <- mis(actual, lower = apply(object$distribution, 2, quantile, alpha/2), upper = apply(object$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
  } else {
    m_mis <- as.numeric(NA)
  }
  m_crps <- crps(actual, object$distribution)
  data.frame("h" = n, "MAPE" = m_mape, "MASE" = m_mase, "MSLRE" = m_mslre, "BIAS" = m_bias, "MIS" = m_mis, "CRPS" = m_crps)
}

#' @method tsmetrics tsets.estimate
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.tsets.estimate = function(object, ...)
{
  # residuals diagnostics
  fx <- fitted(object)
  ac <- xts(object$spec$target$y_orig, object$spec$target$index)
  acfx <- na.omit(cbind(ac, fx))
  actual <- as.numeric(acfx[,1])
  fted <- as.numeric(acfx[,2])
  np <- NROW(object$model$setup$parmatrix[which(object$model$setup$parmatrix[,"required"] == 1),])
  np <- np + 1
  if (!is.null(object$model$setup$theta_type_model) && object$model$setup$theta_type_model) {
    np <- np - 1
  }
  ny <- NROW(acfx)
  if (substr(object$model$setup$model,3,3) != "N") {
    frequency <- object$spec$target$frequency
  } else{
    frequency <- 1
  }
  m_mape <- mape(actual, fted)
  m_mase <- mase(actual, fted, original_series = actual, frequency = frequency)
  m_mslre  <- mslre(actual, fted)
  m_bias <- bias(actual, fted)
  ll   <- object$model$loglik
  lik  <- -0.5 * ll
  aic  <- ll + 2 * np
  bic  <- ll + log(ny) * np
  aicc <- aic + 2 * np * (np + 1) / (ny - np - 1)
  data.frame("n" = ny, "no.pars" = np - 1, "LogLik" = lik, "AIC" = aic, "BIC" = bic, "AICc" = aicc, "MAPE" = m_mape, "MASE" = m_mase,
             "MSLRE" = m_mslre, "BIAS" = m_bias)
}

#' Akaike's An Information Criterion
#'
#' @description Extract the AIC from an estimated model.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param ... not currently used.
#' @param k the penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @return a numeric value.
#' @aliases AIC
#' @method AIC tsets.estimate
#' @rdname AIC
#' @export
#'
#'
AIC.tsets.estimate <- function(object, ..., k = 2)
{
  np <- NROW(object$model$setup$parmatrix[which(object$model$setup$parmatrix[,"required"] == 1),])
  np <- np + 1
  if (!is.null(object$model$setup$theta_type_model) && object$model$setup$theta_type_model) {
    np <- np - 1
  }
  ll   <- object$model$loglik
  lik  <- -0.5 * ll
  aic  <- ll + k * np
  return(aic)
}

#' Extract Model Coefficients
#'
#' @description Extract the estimated coefficients of a model.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param ... not currently used.
#' @return a numeric named vector.
#' @aliases coef
#' @method coef tsets.estimate
#' @rdname coef
#' @export
#'
#'
coef.tsets.estimate = function(object, ...)
{
  v <- object$model$setup$parmatrix[which(object$model$setup$parmatrix[,"required"] == 1),1]
  v <- c(v, "sigma" = object$model$setup$parmatrix["sigma",1])
  return(v)
}

#' Model Estimation Summary
#'
#' @description Summary method for class \dQuote{tsets.estimate}
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... not currently used.
#' @return A printout of the parameter summary, model type and some model metrics.
#' When estimated using autodiff, the standard errors, t-values and p-values will
#' also be printed. In this case, if the parameters are close to their upper or
#' lower bounds then it is very likely that these values will be NaN.
#' @aliases summary
#' @method summary tsets.estimate
#' @rdname summary
#' @export
#'
#'
summary.tsets.estimate = function(object, digits = 4, ...)
{
  parmatrix <- object$model$setup$parmatrix
  pmatrix <- parmatrix
  pmatrix   <- rbind(pmatrix,parmatrix["sigma",,drop = FALSE])
  pmatrix   <- cbind(pmatrix,data.frame("." = rep(NA, nrow(pmatrix))))
  pnames    <- rownames(pmatrix)

  # augment with names
  if (object$model$setup$include_xreg == 1) {
    xreg_names <- paste0("Regressor[",colnames(object$spec$xreg$xreg),"]")
    pmatrix[which(grepl("rho",pnames)),"."] <- xreg_names
  }

  if (object$model$setup$include_trend == 1) {
    if (!is.null(object$model$setup$theta_type_model) && object$model$setup$theta_type_model) {
      trend_names <- "State[Slope-init]"
      pmatrix[which(grepl("b0",pnames)),"."] <- trend_names
    } else {
      trend_names <- c("State[Slope-coef]","State[Slope-init]")
      pmatrix[which(grepl("beta|b0",pnames)),"."] <- trend_names
    }

    if (object$spec$model$include_damped == 1) {
      pmatrix[which(grepl("phi",pnames)),"."] <- "State[Slope-dampening]"
    }
  }

  if (object$model$setup$include_seasonal == 1) {
    pmatrix[which(grepl("gamma",pnames)),"."] <- "State[Seasonal-coef]"
    pmatrix[which(grepl("[s][0-9]",pnames)),"."] <- rep("State[Seasonal-init]",object$spec$seasonal$frequency - 1)
  }

  if (any(pnames == "delta")) {
    pmatrix[which(grepl("delta",pnames)),"."] <- "State[Seasonal-power]"
  }

  if (any(pnames == "theta")) {
    pmatrix[which(grepl("theta",pnames)),"."] <- "State[Trend-power]"
  }

  pmatrix[which(grepl("alpha",pnames)),"."] <- "State[Level-coef]"
  pmatrix[which(grepl("l0",pnames)),"."] <- "State[Level-init]"
  pmatrix[which(grepl("sigma",pnames)),"."] <- "Observation[Standard Deviation]"
  # only using what we estimated
  pmatrix <- pmatrix[which(pmatrix[,"estimate"] == 1), ]
  if (object$spec$model$include_damped == 1) {
    modelx <- paste0(substr(object$spec$model$model,1,1),substr(object$spec$model$model,2,2),"d",substr(object$spec$model$model,3,3))
  } else {
    modelx <- object$spec$model$model
  }

  if (object$spec$model$type == 4) {
    cat("\nPower ETS Model [",modelx,"]")
  } else {
    cat("\nETS Model [",modelx,"]")
  }
  printout <- data.frame("Parameter" = rownames(pmatrix), "Description" = pmatrix[,"."], "Est[Value]" = pmatrix[,1], check.names = FALSE)
  if (!is.null(object$hess)) {
    S <- try(suppressWarnings(.make_standard_errors(object)), silent = TRUE)
    if (!inherits(S,'try-error')) {
      printout <- cbind(printout, S)
    }
  }
  print(kable(printout, right = FALSE, digits = digits, row.names = FALSE, format = "simple"))

  mtrcs <- tsmetrics(object)
  print(kable(data.frame(LogLik = as.numeric(sprintf(mtrcs$LogLik,fmt = "%.4f")),
                         AIC = as.numeric(sprintf(mtrcs$AIC,fmt = "%.2f")),BIC = as.numeric(sprintf(mtrcs$BIC,fmt = "%.2f")),
                   AICc = as.numeric(sprintf(mtrcs$AICc,fmt = "%.2f"))), format = "simple",row.names = FALSE, right = FALSE))
  print(kable(data.frame(MAPE = as.numeric(sprintf(mtrcs$MAPE,fmt = "%.4f")), MASE = as.numeric(sprintf(mtrcs$MASE,fmt = "%.4f")),
                         MSLRE = as.numeric(sprintf(mtrcs$MSLRE,fmt = "%.4f")), BIAS = as.numeric(sprintf(mtrcs$BIAS,fmt = "%.4f"))), format = "simple",row.names = FALSE, right = FALSE))

  if (!is.null(object$selection)) {
    cat("\nModel Ranking\n")
    cat("--------------\n")
    print(kable(object$selection, digits = digits, row.names = FALSE, format = "simple"))
  }
  return(invisible(object))
}

.make_standard_errors <- function(object)
{
  pmatrix <- object$model$setup$parmatrix
  pars <- pmatrix[which(pmatrix[,"estimate"] == 1),1]
  H <- object$hess
  se <- sqrt(diag(solve(H)))
  tvalues <- pars/se
  pvalues <- 2*(1 - pnorm(abs(tvalues)))
  return(data.frame("Std. Error" = se,"t value" = tvalues, "Pr(>|t|)" = pvalues, check.names = FALSE))
}


.tables.tsets.estimate <- function(object, digits = 4, ...)
{
  parmatrix <- object$model$setup$parmatrix
  pmatrix <- parmatrix[which(parmatrix[,"estimate"] == 1),]
  pmatrix <- rbind(pmatrix,parmatrix["sigma",,drop = FALSE])
  pmatrix <- cbind(pmatrix,data.frame("." = rep(NA, nrow(pmatrix))))
  pnames <- rownames(pmatrix)

  # augment with names
  if (object$model$setup$include_xreg == 1) {
    xreg_names <- paste0("Regressor[",colnames(object$spec$xreg$xreg),"]")
    pmatrix[which(grepl("rho",pnames)),"."] <- xreg_names
  }

  if (object$model$setup$include_trend == 1) {
    if (!is.null(object$model$setup$theta_type_model) && object$model$setup$theta_type_model) {
      trend_names <- "State[Slope-init]"
      pmatrix[which(grepl("b0",pnames)),"."] <- trend_names
    } else {
      trend_names <- c("State[Slope-coef]","State[Slope-init]")
      pmatrix[which(grepl("beta|b0",pnames)),"."] <- trend_names
    }

    if (object$spec$model$include_damped == 1) {
      pmatrix[which(grepl("phi",pnames)),"."] <- "State[Slope-dampening]"
    }
  }

  if (object$model$setup$include_seasonal == 1) {
    pmatrix[which(grepl("gamma",pnames)),"."] <- "State[Seasonal-coef]"
    pmatrix[which(grepl("[s][0-9]",pnames)),"."] <- rep("State[Seasonal-init]",object$spec$seasonal$frequency - 1)
  }

  if (any(pnames == "delta")) {
    pmatrix[which(grepl("delta",pnames)),"."] <- "State[Seasonal-power]"
  }

  if (any(pnames == "theta")) {
    pmatrix[which(grepl("theta",pnames)),"."] <- "State[Trend-power]"
  }

  pmatrix[which(grepl("alpha",pnames)),"."] <- "State[Level-coef]"
  pmatrix[which(grepl("l0",pnames)),"."] <- "State[Level-init]"
  pmatrix[which(grepl("sigma",pnames)),"."] <- "Observation[Standard Deviation]"

  if (object$spec$model$include_damped == 1) {
    modelx <- paste0(substr(object$spec$model$model,1,1),substr(object$spec$model$model,2,2),"d",substr(object$spec$model$model,3,3))
  } else {
    modelx <- object$spec$model$model
  }

  if (object$spec$model$type == 4) {
    modelx <- paste0("Power ETS Model [",modelx,"]")
  } else {
    modelx <- paste0("ETS Model [",modelx,"]")
  }

  if (any(grepl("[rho][0-9]",rownames(pmatrix)))) {
    pmatrixrho <- pmatrix[which(grepl("[rho][0-9]",rownames(pmatrix))),]
    pmatrixnonrho <- pmatrix[-which(grepl("[rho][0-9]",rownames(pmatrix))),]
    p1 <- data.frame("Parameter" = rownames(pmatrixnonrho),"Description" = pmatrixnonrho[,"."], "Est[Value]" = pmatrixnonrho[,1], check.names = FALSE)
    p2 <- data.frame("Parameter" = rownames(pmatrixrho),"Description" = pmatrixrho[,"."], "Est[Value]" = pmatrixrho[,1], check.names = FALSE)
    printout <- list(p1, p2)
    n <- 2
  } else {
    printout <- data.frame("Parameter" = rownames(pmatrix),"Description" = pmatrix[,"."], "Est[Value]" = pmatrix[,1], check.names = FALSE)
    n <- 1
  }
  m <- tsmetrics(object)
  infomatrix <- data.frame(AIC = as.numeric(sprintf(m$AIC,fmt = "%.2f")),BIC = as.numeric(sprintf(m$BIC,fmt = "%.2f")), AICc = as.numeric(sprintf(m$AICc,fmt = "%.2f")))
  mmatrix <- data.frame(MAPE = as.numeric(sprintf(m$MAPE,fmt = "%.4f")), MASE = as.numeric(sprintf(m$MASE,fmt = "%.4f")),
                       MSLRE = as.numeric(sprintf(m$MSLRE,fmt = "%.4f")),BIAS = as.numeric(sprintf(m$BIAS,fmt = "%.4f")))
  if (!is.null(object$selection)) {
    selection <- object$selection
  } else {
    selection <- NULL
  }

  return(list(model = modelx, params = printout, n = n, info = infomatrix, metrics = mmatrix, selection = selection))
}

#' Estimation Summary Report (pdf format)
#'
#' @description Generates a pdf summary of the estimated model.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param output_dir directory where the output is written and the object saved.
#' @param args additional list of arguments used in the generation of the report.
#' Only name is currently used to display the name of the data series.
#' @param ... not currently used.
#' @return A pdf file.
#' @aliases tsreport
#' @method tsreport tsets.estimate
#' @rdname tsreport
#' @export
#'
#'
tsreport.tsets.estimate <- function(object, output_dir = "/", args = list(name = NULL), ...)
{
  saveRDS(object, file = paste0(output_dir,"/tsets_estimate_tmp.rds"))
  fname <- args$name
  if (is.null(fname)) fname <- paste0(fname,"_tsets_estimation_report_",format(as.Date(Sys.Date()),"%Y%m%d")) else fname <- paste0("tsets_estimation_report_",format(as.Date(Sys.Date()),"%Y%m%d"))
  suppressWarnings(rmarkdown::render(file.path(find.package('tsets'),'reports/estimation_report.Rmd'),
                                     output_dir = output_dir,
                                     output_file = paste0(fname,".pdf"),
                                     intermediates_dir = tempdir(),
                                     output_format = pdf_document2(toc = FALSE),
                                     params = list(dir = as.character(output_dir), name = args$name), quiet = TRUE))
}


#' Model Decomposition
#'
#' @description Decomposes the estimated model or prediction into its component
#' parts (states).
#' @param object an object of class \dQuote{tsets.estimate} or \dQuote{tsets.predict}
#' @param simplify whether to return the components as Trend (Level + Slope), Seasonal,
#' X and Irregular.
#' @param ... not currently used.
#' @return For the estimated object, returns an xts matrix of the state components
#' (including error). For the predicted object, a list with the simulated object
#' which includes the predictive distribution matrix of the state components and
#' their estimated values, inheriting class \dQuote{tsmodel.predict} (for pretty
#' plotting functionality).
#' @aliases tsdecompose
#' @method tsdecompose tsets.estimate
#' @rdname tsdecompose
#' @export
#'
#'
tsdecompose.tsets.estimate <- function(object, simplify = FALSE, ...)
{
  d <- object$spec$target$index
  f <- fitted(object)
  # level and slope
  idx <- 1:(nrow(object$model$states) - 1)
  Level <- xts(object$model$states[idx,"Level"],d)
  colnames(Level) <- "Level"
  if (simplify) {
    Trend <- Level
  }
  if (object$spec$model$include_trend == 1) {
    Slope <- xts(object$model$states[idx,"Trend"],d)
    colnames(Slope) <- "Slope"
    if (simplify) {
      Trend <- Trend + Slope
    }
  } else {
    Slope <- NULL
  }
  # seasonal
  if (object$spec$model$include_seasonal == 1) {
    Seasonal <- xts(object$model$states[idx, paste0("S",object$spec$seasonal$frequency - 1)],d)
    colnames(Seasonal) <- "Seasonal"
  } else {
    Seasonal <- NULL
  }
  # external
  if (object$spec$model$include_xreg == 1) {
    k     <- ncol(object$spec$xreg$xreg)
    xreg  <- object$spec$xreg$xreg
    rho   <- matrix(object$model$setup$parmatrix[paste0("rho",1:k),1], ncol = 1)
    X     <- xts(xreg %*% rho,d)
    colnames(X) <- "X"
  } else {
    X <- NULL
  }
  # residuals
  Irregular <- xts(object$model$residuals, d)
  colnames(Irregular) <- "Irregular"
  if (simplify) {
    return(list(Trend = Trend, Seasonal = Seasonal, X = X, Irregular = Irregular))
  } else {
    return(list(Level = Level, Slope = Slope, Seasonal = Seasonal, X = X, Irregular = Irregular))
  }
}

#' @method tsdecompose tsets.predict
#' @rdname tsdecompose
#' @export
#'
tsdecompose.tsets.predict <- function(object, simplify = FALSE, ...)
{
  if (simplify) {
    d <- object$decomposition
    cnames <- names(d)
    Trend <- d$Level
    if (any(cnames %in% "Slope")) {
      Trend$distribution <- Trend$distribution + d$Slope$distribution
      Trend$original_series <- Trend$original_series + d$Slope$original_series
    }
    if (any(cnames %in% "Seasonal")) {
      Seasonal <- d$Seasonal
    } else {
      Seasonal <- NULL
    }
    if (any(cnames %in% "X")) {
      X <- d$X
    } else {
      X <- NULL
    }
    Irregular <- d$Error
    d <- list(Trend = Trend, Seasonal = Seasonal, X = X, Irregular = Irregular)
  } else {
    d <- object$decomposition
    return(d)
  }
}

#' Model Specification Extractor
#'
#' @description Extracts a model specification (class \dQuote{tsets.spec}) from
#' an object of class \dQuote{tsets.estimate}.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param y an optional new xts vector.
#' @param lambda an optional lambda parameter for the Box Cox transformation (if
#' previously used).
#' @param xreg an optional matrix of regressors.
#' @param ... not currently used.
#' @note This function is used by other functions in the package such as the
#' backtest which requires rebuilding the specification for each re-estimation
#' step with updated data but keeping all else equal.
#' @return An object of class \dQuote{tsets.spec}.
#' @aliases tsspec
#' @method tsspec tsets.estimate
#' @rdname tsspec
#' @export
#'
#'
tsspec.tsets.estimate <- function(object, y = NULL, lambda = NULL, xreg = NULL, ...)
{
    if (is.null(y)) {
        y <- xts(object$spec$target$y.orig, object$spec$target$index)
    }
    if (!is.null(xreg)) {
        xreg <- coredata(xreg)
        if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
    } else {
        if (object$spec$xreg$include_xreg) {
        xreg <- object$spec$xreg$xreg
        if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
        } else {
            xreg <- NULL
        }
    }
    if (!is.null(object$spec$transform)) {
        transformation <- object$spec$transform$name
        lower <- object$spec$transform$lower
        upper <- object$spec$transform$upper
        if (is.null(lambda)) {
            lambda <- object$spec$transform$lambda
        }
    } else {
        transformation <- NULL
        lower <- NULL
        upper <- NULL
        lambda <- NULL
    }


  spec <- ets_modelspec(y, model = object$spec$model$model, damped = object$spec$model$damped,
                        power = object$spec$model$power, xreg = xreg,
                        frequency = object$spec$target$frequency, transformation = transformation,
                        lambda = lambda, lower = lower, upper = upper,
                        normalized_seasonality = object$spec$model$normalized_seasonality,
                        scale = attr(object$spec$target$scaler, "scale"),
                        seasonal_init = object$spec$model$seasonal_init)
  return(spec)
}
