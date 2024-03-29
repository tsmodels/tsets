#' Walk Forward Model Backtest
#'
#' @description Generates an expanding window walk forward backtest.
#' @param object an object of class \dQuote{tsets.spec}.
#' @param start numeric data index from which to start the backtest.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param estimate_every number of periods at which the model is re-estimated
#' and new predictions are generated (defaults to 1).
#' @param FUN optional function which is applied across all horizons for each
#' draw (i.e. operating on each row of the distribution which represents a
#' single path). For example, using the max function will return the distribution
#' of the maximum across all horizons, whereas sum (for flow variables) would
#' represent the aggregate value distribution. The P50 of this distribution is
#' returned and aligned with the terminal horizon for each re-estimation period,
#' and if alpha is not NULL, then the quantiles of this distributions with
#' respect to the coverage (alpha) chosen.
#' @param alpha optional numeric vector of coverage rates for which to calculate
#' the quantiles.
#' @param solver solver to use.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsetsad package.
#' @param autoclean whether to perform automatic cleaning on the training data
#' prior to prediction as per the \sQuote{auto_clean} function in the tsaux package.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... additional arguments passed to the \dQuote{auto_clean} function.
#' @return A list with the following data.tables:
#' \itemize{
#' \item prediction : the backtest table with forecasts and actuals
#' \item metrics: a summary performance table showing metrics by
#' forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
#' }
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package.
#' @aliases tsbacktest
#' @method tsbacktest tsets.spec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.tsets.spec <- function(object, start = floor(length(object$target$y_orig))/2, end = length(object$target$y_orig),
                                  h = 1, estimate_every = 1, FUN = NULL, alpha = NULL, solver = "nlminb",
                                  autodiff = FALSE, autoclean = FALSE, trace = FALSE, ...)
{
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
    if (estimate_every != 1) {
        estimate_every <- max(1, as.integer(estimate_every))
        ns <- length(seqdates)
        seqdates <- seqdates[seq(1, ns, by = estimate_every)]
    }
    if (!is.null(alpha)) {
        if (any(alpha <= 0)) {
            stop("\nalpha must be strictly positive")
        }
        if (any(alpha >= 1)) {
            stop("\nalpha must be less than 1")
        }
        quantiles <- as.vector(sapply(1:length(alpha), function(k) c(alpha[k]/2, 1 - alpha[k]/2)))
    } else {
        quantiles <- NULL
    }
    # setup backtest indices
    horizon <- sapply(1:length(seqdates), function(i){
        min(h, elapsed_time(index(data), index(data)[end], seqdates[i]))
    })
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    extra_args <- list(...)
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
        p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(ytest))
        if (!is.null(quantiles)) {
            qp <- apply(p$distribution, 2, quantile, quantiles)
            if (length(quantiles) == 1) {
                qp <- matrix(qp, ncol = 1)
            } else{
                qp <- t(qp)
            }
            colnames(qp) <- paste0("P", round(quantiles*100))
        }
        out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]),
                          "horizon" = 1:horizon[i],
                          "size" = rep(nrow(ytrain), horizon[i]),
                          "forecast_dates" = as.character(index(ytest)),
                          "forecast" = as.numeric(p$mean), "actual" = as.numeric(ytest))
        if (!is.null(quantiles)) out <- cbind(out, qp)
        if (!is.null(FUN)) {
            pd <- apply(p$distribution, 1, FUN)
            funp <- data.table(estimation_date = seqdates[i], horizon = horizon[i], fun_P50 = median(pd), fun_actual = FUN(as.numeric(ytest)))
            if (!is.null(quantiles)) {
                qp <- matrix(quantile(pd, quantiles), nrow = 1)
                colnames(qp) <- paste0("fun_P", round(quantiles*100,1))
                funp <- cbind(funp, qp)
            }
            out <- merge(out, funp, by = c("estimation_date","horizon"), all = TRUE)
        }
        return(out)
    }, future.packages = c("tsmethods","tsaux","xts","tsets","data.table"), future.seed = TRUE)
    b <- eval(b)
    b <- rbindlist(b)
    data_name <- "y"
    actual <- NULL
    forecast <- NULL
    metrics <- b[,list(variable = data_name, MAPE = mape(actual, forecast), MSLRE = mslre(actual, forecast),
                       BIAS = bias(actual, forecast),
                       n = .N), by = "horizon"]
    if (!is.null(alpha)) {
        q_names <- matrix(paste0("P", round(quantiles*100)), ncol = 2, byrow = TRUE)
        q <- do.call(cbind, lapply(1:length(alpha), function(i){
            b[,list(mis = mis(actual, get(q_names[i,1]), get(q_names[i,2]), alpha[i])), by = "horizon"]
        }))
        q <- q[,which(grepl("mis",colnames(q))), with = FALSE]
        colnames(q) <- paste0("MIS[",alpha,"]")
        metrics <- cbind(metrics, q)
    }
    return(list(prediction = b, metrics = metrics))
}


