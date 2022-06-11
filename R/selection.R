#' Automatic Model Selection
#'
#' @description Automatic model selection based on choice of criteria.
#' @details The function estimates a number of different models from the class
#' of valid ETS models in the packages, with and without dampening, with and
#' without a power term (for MAM and MAN models).
#' @param y an xts vector.
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
#' @param metric the selection metric to use. Valid metrics are \sQuote{AIC},
#' \sQuote{BIC}, \sQuote{AICc}, \sQuote{MASE} and \sQuote{MAPE}. If lambda is
#' not NULL, then those models which admit a Box-Cox transformation (additive
#' error models) will not be comparable with the other models.
#' @param additive_only whether to limit to additive models only.
#' @param solver the solver to use for estimation.
#' @param control the solver control parameters.
#' @param power_model whether to include the power MAM models.
#' @param include_damped whether to include damped models in the selection.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param return_table whether to return the table with the enumerated options,
#' ranked by metric,for each combination of those options used.
#' @param scale whether to rescale the data using y/max(y) (only for additive models).
#' This sometimes helps in the optimization.
#' @param seasonal_init whether the initial seasonal states are estimated or
#' fixed (set to a backcast approximation).
#' @param autodiff whether to use automatic differentiation
#' (see \code{\link{estimate.tsets.spec}}).
#' @param ... not used.
#' @return An object of class \dQuote{tsets.estimate} which also inherits
#' class \dQuote{tsets.select}
#' @note The function can use parallel functionality as long as the user has set up a
#' \code{\link[future]{plan}} using the future package.
#' @aliases auto_ets
#' @rdname auto_ets
#' @export
#'
auto_ets = function(y, xreg = NULL, transformation = NULL, lambda = NULL, lower = 0, upper = 1,
                    metric = "AIC", frequency = NULL, normalized_seasonality = TRUE,
                    additive_only = FALSE, solver = "nlminb",
                    control = list(trace = 0, maxit = 1000), power_model = FALSE,
                    include_damped = TRUE, trace = FALSE, return_table = FALSE,
                    scale = FALSE, seasonal_init = "fixed", autodiff = TRUE, ...)
{
  # sanity check of model input
  valid_criteria <- c("AIC","BIC","AICc","MAPE","MASE","MSLRE")
  if (!(metric %in% valid_criteria)) {
    stop("\nvalid criteria are AIC, BIC, AICc, MAPE, MASE and MSLRE")
  }
  if ((metric == "AIC" | metric == "BIC" | metric == "AICc") && !additive_only) {
    if (!is.null(transformation)) warning("\ninvalid selection metric when transformation is not NULL (models will not be comparable)")
  }
  if (any(y < 0, na.rm = T) && !additive_only) {
    warning("the data contains negative values; not suitable for multiplicative models!")
  }
  valid_models <- c("AAA","AAN","ANN","ANA","Theta")
  if (!additive_only) {
    valid_models <- c(valid_models,"MMM","MMN","MNM","MNN","MAM","MAN")
    # add the power models
    valid_models <- c(valid_models,"MAM","MAN")
  }
  n_valid_models <- length(valid_models)

  sgrid <- data.frame(model = valid_models, power = 0, stringsAsFactors = FALSE)

  # process additive and multiplicative models separately
  if (!additive_only & power_model) {
    sgrid[which(sgrid$model == "MAM")[2],"power"] <- 1
    sgrid[which(sgrid$model == "MAN")[2],"power"] <- 1
  }

  sgrid <- rbind( cbind(sgrid, data.frame(damped = rep(0, nrow(sgrid)))), cbind(sgrid, data.frame(damped = rep(1,nrow(sgrid)))) )
  sgrid[which(substr(sgrid$model,2,2) == "N"), "damped"] <- 0
  if (!include_damped) {
    sgrid <- sgrid[-which(sgrid$damped == 1),]
  }
  ngrid <- NROW(sgrid)

  unid <- sapply(1:ngrid, function(i) paste0(sgrid$model[i], sgrid$power[i], sgrid$damped[i]))
  exc <- which(duplicated(unid))
  sgrid <- sgrid[-exc,]

  rownames(sgrid) <- NULL
  ngrid <- NROW(sgrid)

  # reduce models if frequency == 1
  if (is.null(frequency) || frequency == 1) {
    exc <- which(substr(sgrid$model,3,3) != "N" & sgrid$model != "Theta")
    sgrid <- sgrid[-exc,]
    ngrid <- NROW(sgrid)
    warnings("\nnot estimating seasonal models (frequency = 1 or NULL)")
  }
  if (trace) {
    prog_trace <- progressor(ngrid)
  }
  # ToDo: retain only coefs and metrics. Rerun for optimal model once more at end  (similar to tsissm seelection)
  models %<-% future_lapply(1:ngrid, function(i) {
      if (trace) prog_trace()
      if (substr(sgrid[i,'model'],1,1) == "M") {
          lambda <- NULL
          trm <- NULL
      } else {
          trm <- transformation[1]
      }
      damped_option_i <- as.logical(sgrid[i,"damped"])
      power_option_i <- as.logical(sgrid[i,"power"])
      spec <- suppressWarnings(ets_modelspec(y, model = sgrid[i,'model'], transformation = trm,
                                             lambda =  lambda, damped = damped_option_i, power = power_option_i,
                                             xreg = xreg, frequency = frequency, normalized_seasonality = normalized_seasonality,
                                             scale = scale, seasonal_init = seasonal_init, lower = lower, upper = upper))
      mod <- estimate(spec, solver = solver, control = control, autodiff = autodiff)
      return(mod)
  }, future.packages = c("tsets","xts","tsaux"))
  models <- eval(models)
  converge <- data.frame("Converged" = as.logical(1 - sapply(models, function(x) x$opt$convergence)), stringsAsFactors = FALSE)
  metricsmat <- do.call(rbind, lapply(1:length(models), function(i) tsmetrics(models[[i]])))
  parmat <- do.call.fast(rbind, lapply(1:length(models), function(i){
    models[[i]]$model$setup$parmatrix[c("alpha","beta","gamma","phi","theta","delta","sigma"),1]
  }))
  outmat <- cbind(sgrid, parmat, metricsmat, converge)
  crt <- sort.int(outmat[,metric], index.return = TRUE)
  outmat <- outmat[crt$ix,]
  rownames(outmat) <- NULL
  optmodel <- models[[crt$ix[1]]]
  if (return_table) {
    optmodel$selection <- outmat
  }
  class(optmodel) <- c(class(optmodel), "tsets.select")
  return(optmodel)
}
