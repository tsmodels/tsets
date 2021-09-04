tsbacktest.tsets.spec <- function(object, start = floor(length(object$target$y_orig))/2, end = length(object$target$y_orig),
                                  h = 1, estimate_every = 1, FUN = NULL, alpha = NULL, cores = 1, data_name = "y", save_output = FALSE, 
                                  save_dir = "~/tmp/", solver = "nlminb", autodiff = FALSE, trace = FALSE, ...)
{
    if (save_output) {
        if (is.null(save_dir)) {
            stop("save_dir cannot be NULL when save.output is TRUE")
        }
        if (!dir.exists(save_dir)) {
            stop("save_dir does not exit. Create first and then resubmit")
        }
    }
    data <- xts(object$target$y_orig, object$target$index)
    model <- object$model$model
    damped <- object$model$damped
    if (!is.null(object$model$power)) {
        power <- object$model$power
    } else{
        power <- FALSE
    }
    transform <- object$transform
    lambda <- transform$lambda
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
    i <- 1
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    clusterExport(cl, "FUN", envir = environment())
    if (trace) {
        iterations <- length(seqdates)
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    b <- foreach(i = 1:length(seqdates), .packages = c("tsmethods","tsaux","xts","tsets","data.table"), .options.snow = opts, .combine = rbind) %dopar% {
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
        spec <- ets_modelspec(ytrain, model = model, damped = damped, power = power, xreg = xreg_train, 
                              frequency = frequency, transform = transform, 
                              normalized_seasonality = normalized_seasonality, seasonal_init = seasonal_init)
        mod <- estimate(spec, solver = solver, autodiff = autodiff)
        p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(ytest))
        if (save_output) {
            saveRDS(mod, file = paste0(save_dir,"/model_", seqdates[i], ".rds"))
            saveRDS(p, file = paste0(save_dir,"/predict_", seqdates[i], ".rds"))
        }
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
    }
    stopCluster(cl)
    if (trace) {
        close(pb)
    }
    if (is.null(data_name)) data_name <- "y"
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


