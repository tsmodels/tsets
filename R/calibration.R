tscalibrate.tsets.spec <- function(object, start = floor(length(object$target$y_orig))/2, end = length(object$target$y_orig),
                                   h = 1, nsim = 5000, cores = 1, solver = "nlminb", autodiff = FALSE,
                                   autoclean = FALSE, trace = FALSE, ...)
{
    if (object$model$type == 4) stop("\nno prediction calibration currently available for the power MAM model.")
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
    i <- 1
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (trace) {
        iterations <- length(seqdates)
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    extra_args <- list(...)
    b <- foreach(i = 1:length(seqdates), .packages = c("tsmethods","tsaux","xts","tsets","data.table","bootstrap"), .options.snow = opts, .combine = rbind) %dopar% {
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
            sigma_d <- sqrt(analytic_moments(mod, h = horizon[i])$var)
        } else {
            er <- as.numeric(ytest) - as.numeric(p$mean)
            sigma_d <- sqrt(analytic_moments(mod, h = horizon[i])$var)
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
    }
    stopCluster(cl)
    if (trace) {
        close(pb)
    }
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
