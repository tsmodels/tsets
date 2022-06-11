#' Walk Forward Calibration of Distribution Variance (Experimental)
#'
#' @description Using an expanding window walk forward backtest, provides a
#' calibrated set of adjustment values which can be used in the predict function
#' for scaling the multi-horizon distributional standard deviation.
#' @param object an object of class \dQuote{tsets.spec}.
#' @param start numeric data index from which to start the backtest calibration.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param nsim number of samples to draw for the simulated prediction distribution.
#' @param solver solver to use.
#' @param autodiff whether to use automatic differentiation for estimation. This
#' makes use of the tsetsad package.
#' @param autoclean whether to perform automatic cleaning on the training data
#' prior to prediction as per the \sQuote{auto_clean} function in the tsaux package.
#' @param trace whether to show the progress bar. The user is expected to have set
#' up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... Additional arguments passed to the \dQuote{auto_clean} function.
#' @return A list with the following data.tables:
#' \itemize{
#' \item sample : multi-horizon white noise sample
#' \item sigma_scale: The scaling values per horizon which can be used to scale
#' the predictive distribution standard deviation.
#' \item sd_sigma_scale: The standard deviation of the sigma scaling factor
#' }
#' @note The walk forward predictive errors per horizon are used to estimate the
#' realized variance against the expected variance (sigma_scale). A jackknife
#' procedure is used for estimating the mean of this value. The output also
#' includes samples from a kernel fitted on the whitened errors per horizon
#' which can then be used as inputs to the prediction function. The function can
#' use parallel functionality as long as the user has set up a
#' \code{\link[future]{plan}} using the future package.
#' @aliases tscalibrate
#' @method tscalibrate tsets.spec
#' @rdname tscalibrate
#' @export
#'
tscalibrate.tsets.spec <- function(object, start = floor(length(object$target$y_orig))/2, end = length(object$target$y_orig),
                                   h = 1, nsim = 5000, solver = "nlminb", autodiff = TRUE,
                                   autoclean = FALSE, trace = FALSE, ...)
{
    if (object$model$type == 4) stop("\nno prediction calibration currently available for the power MAM model.")
    if (object$model$class > 3) stop("\nno prediction calibration currently available for models of class 4 (MMM, MMN).")
    data <- xts(object$target$y_orig, object$target$index)
    model <- object$model$model
    damped <- object$model$damped
    if (!is.null(object$model$power)) {
        power <- object$model$power
    } else{
        power <- FALSE
    }
    transform <- object$transform
    if (!is.null(transform)) {
        if (transform$name == "box-cox") {
            if (transform$estimated) {
                lambda <- NA
            } else {
                lambda <- transform$lambda
            }
        } else {
            lambda <- NULL
        }
    } else {
        transform <- list()
        lambda <- NULL
    }
    normalized_seasonality <- object$model$normalized_seasonality
    frequency <- object$target$frequency
    seasonal_init <- object$model$seasonal_init
    if (object$model$include_xreg == 1) {
        use_xreg <- TRUE
        xreg <- xts(object$xreg$xreg, object$target$index)
    } else {
        use_xreg <- FALSE
        xreg <- NULL
    }
    start_date <- index(data)[start]
    end <- min(NROW(data), end)
    end_date <- index(data)[end - 1]
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    elapsed_time <- function(idx, end_date, start_date) {
        min(h, which(end_date == idx) - which(start_date == idx))
    }
    # setup backtest indices
    horizon <- sapply(1:length(seqdates), function(i){
        min(h, elapsed_time(index(data), index(data)[end], seqdates[i]))
    })
    extra_args <- list(...)
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    b %<-% future_lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        ytrain <- data[paste0("/", seqdates[i])]
        ix <- which(index(data) == seqdates[i])
        ytest <- data[(ix + 1):(ix + horizon[i])]
        if (use_xreg) {
            xreg_train <- xreg[index(ytrain)]
            xreg_test <- xreg[index(ytest)]
        } else {
            xreg_train <- NULL
            xreg_test <- NULL
        }
        if (autoclean) {
            args_x <- c(list(y = ytrain), list(frequency = frequency), list(lambda = lambda), extra_args)
            ytrain <- do.call(auto_clean, args = args_x, quote = TRUE)
        }
        spec <- ets_modelspec(ytrain, model = model, damped = damped, power = power, xreg = xreg_train,
                              frequency = frequency, lambda = lambda, transformation = transform$name,
                              lower = transform$lower, upper = transform$upper,
                              normalized_seasonality = normalized_seasonality, seasonal_init = seasonal_init)
        mod <- estimate(spec, solver = solver, autodiff = autodiff)
        p <- predict(mod, h = horizon[i], newxreg = xreg_test, nsim = nsim, forc_dates = index(ytest))
        quantiles <- seq(0.1, 1, by = 0.1)
        qp <- apply(p$distribution, 2, quantile, quantiles)
        qp <- t(qp)
        colnames(qp) <- paste0("P", round(quantiles*100))
        if (!is.null(transform$name)) {
            dist <- do.call(cbind, lapply(1:ncol(p$distribution), function(j) mod$spec$transform$transform(p$distribution[,j], lambda = mod$spec$transform$lambda)))
            ac <- mod$spec$transform$transform(ytest, lambda = mod$spec$transform$lambda)
            er <- unname(as.numeric(ac - colMeans(dist)))
            sigma_d <- sqrt(tsmoments(mod, h = horizon[i], newxreg = xreg_test)$var)
        } else {
            er <- as.numeric(ytest) - as.numeric(p$mean)
            sigma_d <- sqrt(tsmoments(mod, h = horizon[i], newxreg = xreg_test)$var)
        }
        # output the calibration table for the binned quantiles
        out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]),
                          "horizon" = 1:horizon[i],
                          "size" = rep(nrow(ytrain), horizon[i]),
                          "forecast_dates" = as.character(index(ytest)),
                          "forecast" = as.numeric(p$mean),
                          "actual" = as.numeric(ytest),
                          "error"  = er,
                          "sigma" = sigma_d)
        out <- cbind(out, qp)
        return(out)
    }, future.packages = c("tsmethods","tsaux","xts","tsets","data.table","bootstrap"), future.seed = TRUE)
    b <- eval(b)
    b <- rbindlist(b)
    error <- NULL
    bx <- b[,list(scaling = abs(error)/sigma, error = error), by = c("estimation_date","horizon")]
    bs <- b[,list(rel_scale = sigma[1]/sigma, horizon = horizon), by = c("estimation_date")]
    bx <- merge(bx, bs, by = c("estimation_date","horizon"))
    bxe <- dcast(data = bx, estimation_date~horizon, value.var = "error")
    bxe <- as.matrix(bxe[,-1])
    bxe <- scale(bxe)
    if (NROW(na.omit(bxe)) < h) {
        V <- cov(bxe, use = "pairwise.complete")
    } else {
        V <- cov(na.omit(bxe))
    }
    V <- make.positive.definite(V)
    E <- eigen(V)
    W <- diag(1/sqrt(E$value)) %*% t(E$vectors)
    Z <- tcrossprod(bxe, W)
    Z <- sweep(Z, 2, colMeans(Z, na.rm = TRUE), "-")
    sigma_scale <- sapply(1:h, function(i) tail(jackknife(bx[horizon == i]$scaling, median)$jack.value, 1))
    sigma_scale_sigma <- sapply(1:h, function(i) tail(jackknife(bx[horizon == i]$scaling, sd)$jack.value, 1))
    samples <- do.call(cbind, lapply(1:h, function(i){
        kmod <- kde(as.numeric(na.omit(Z[,i])))
        s <- rkde(nsim, kmod)
        return(matrix(s, ncol = 1))
    }))
    return(list(samples = samples, sigma_scale = sigma_scale, sd_sigma_scale = sigma_scale_sigma))
}
