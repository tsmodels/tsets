# automatic model selection
auto_ets = function(y, xreg = NULL, lambda = NA, metric = "AIC", frequency = NULL,
                   normalized_seasonality = TRUE, additive_only = FALSE, cores = NULL, 
                   solver = "nlminb", control = list(trace = 0, maxit = 1000), 
                   power_model = FALSE, verbose = FALSE, retain = 1, scale = FALSE, seasonal_init = "fixed", ...)
{
  # sanity check of model input
  valid_criteria <- c("AIC","BIC","AICc","MAPE","MASE","MSLRE")
  if (!(metric %in% valid_criteria)) {
    stop("\nvalid criteria are AIC, BIC, AICc, MAPE, MASE and MSLRE")
  }
  if (metric == "AIC" | metric == "BIC" | metric == "AICc") {
    if (!is.null(lambda)) warning("\ninvalid selection metric when transform is not NULL (models will not be comparable)")
  }

  if (any(y < 0) && !additive_only) {
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

  # process additive and multiplicative models seperately
  if (!is.null(lambda)) {
    transform <- box_cox(lambda = NA)
    y_t <- transform$transform(y, lambda = NA, frequency = frequency)
    transform$lambda <- attr(y_t,"lambda")
  } else {
    transform <- NULL
  }
  if (!additive_only | power_model) {
    sgrid[which(sgrid$model == "MAM")[2],"power"] <- 1
    sgrid[which(sgrid$model == "MAN")[2],"power"] <- 1
  }

  sgrid <- rbind( cbind(sgrid, data.frame(damped = rep(0, nrow(sgrid)))), cbind(sgrid, data.frame(damped = rep(1,nrow(sgrid)))) )
  sgrid[which(substr(sgrid$model,2,2) == "N"), "damped"] <- 0
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
  i <- 1
  if (is.null(cores)) {
    registerDoSEQ()
    v <- foreach(i = 1:ngrid, .packages = c("tsets","xts","tsaux")) %do% {
      # check if not additive
      if (substr(sgrid[i,'model'],1,1) == "M") {
        lambda <- NULL
      } else {
        lambda <- transform$lambda
      }
      damped_option_i <- as.logical(sgrid[i,"damped"])
      power_option_i <- as.logical(sgrid[i,"power"])

      spec <- suppressWarnings(ets_modelspec(y, model = sgrid[i,'model'], lambda = lambda, damped = damped_option_i, power = power_option_i,
                                            xreg = xreg, frequency = frequency, normalized_seasonality = normalized_seasonality,
                                            scale = scale, seasonal.init = seasonal_init))
      mod <- estimate(spec, solver = solver, control = control)
      return(mod)
    }
  } else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (verbose) {
      pb <- txtProgressBar(max = ngrid, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }
    v <- foreach(i = 1:ngrid, .packages = c("tsets","xts","tsaux"), .options.snow = opts) %dopar% {
      if (substr(sgrid[i,'model'],1,1) == "M") {
        lambda <- NULL
      } else {
        lambda <- transform$lambda
      }
      damped_option_i <- as.logical(sgrid[i,"damped"])
      power_option_i <- as.logical(sgrid[i,"power"])

      spec <- suppressWarnings(ets_modelspec(y, model = sgrid[i,'model'], lambda =  lambda,
                                            damped = damped_option_i, power = power_option_i,
                                            xreg = xreg, frequency = frequency, normalized_seasonality = normalized_seasonality,
                                            scale = scale, seasonal_init = seasonal_init))
      mod <- estimate(spec, solver = solver, control = control)
      return(mod)
    }

    if (verbose) close(pb)

    stopCluster(cl)
  }

  converge <- data.frame("Converged" = as.logical(1 - sapply(v, function(x) x$opt$convergence)), stringsAsFactors = FALSE)
  metricsmat <- do.call(rbind, lapply(1:length(v), function(i) tsmetrics(v[[i]])))
  parmat <- do.call.fast(rbind, lapply(1:length(v), function(i){
    v[[i]]$model$setup$parmatrix[c("alpha","beta","gamma","phi","theta","delta","sigma"),1]
  }))

  outmat <- cbind(sgrid, parmat, metricsmat, converge)

  crt <- sort.int(outmat[,metric], index.return = TRUE)
  outmat <- outmat[crt$ix,]

  rownames(outmat) <- NULL
  optmodel <- v[[crt$ix[1]]]
  optmodel$selection <- outmat

  if (retain > 1) {
    optmodel$retained <- v[crt$ix[2:retain]]
  }
  class(optmodel) <- c(class(optmodel), "tsets.select")
  return(optmodel)
}
