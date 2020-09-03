# nstep ahead in-sample forecast error needed for calculating the h-step covariance matrix
hresiduals.tsets.estimate <- function(object, h = 12, cores = 1, seed = NULL, trace = FALSE, raw = TRUE, ...)
{
    n <- length(object$spec$target$y)
    hfun <- switch(object$spec$model$type,
           "1" = tsets:::hresiduals_aaa_cpp,
           "2" = tsets:::hresiduals_mmm_cpp,
           "3" = tsets:::hresiduals_mam_cpp,
           "4" = tsets:::hresiduals_powermam_cpp)
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (trace) {
        iterations <- n
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    i <- 0
    h_residuals <- foreach(i = 1:(n - 1), .packages = c("tsmethods","tsets","xts","data.table"), .options.snow = opts) %dopar% {
        tmp <- hfun(object, h = h, nsim = 2000, start = i, seed = seed, raw = raw)
        return(tmp)
    }
    if (trace == 1) {
        close(pb)
    }
    stopCluster(cl)
    h_residuals <- rbindlist(h_residuals)
    h_residuals <- dcast(h_residuals, date~horizon, value.var = "error")
    h_residuals <- na.omit(h_residuals)
    return(h_residuals)
}



hresiduals_aaa_cpp <- function(object, h = 12, nsim = 1000, start = 1, seed = NULL, raw = TRUE, ...)
{
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    # truncate h to the available time periods
    coefficient <- object$model$setup$parmatrix[,1]
    frequency <- object$spec$seasonal$frequency
    # last state
    N <- length(object$spec$target$y)
    if ((start + h) > N) {
        h <- N - start
    }
    model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim)
    if (raw) {
        actual <- object$spec$target$y[(start):(start + h - 1)]
    } else {
        actual <- object$spec$target$y_orig[(start):(start + h - 1)]
    }
    # states start one period before (because of initialization seed)
    stateinit <- object$model$states[start,,drop = FALSE]
    pars <- rep(0, 6)
    pars[1] <- stateinit[1,"Level"]
    if (model[1] == 1) pars[2] <- stateinit[1,"Trend"]
    pars[3] <- coefficient["alpha"]
    pars[4] <- coefficient["beta"]
    pars[5] <- coefficient["gamma"]
    pars[6] <- coefficient["phi"]
    pars <- unname(as.numeric(pars))
    if (model[2] == 1) {
        s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
    } else {
        s0 <- rep(0, frequency)
    }
    s0 <- unname(as.numeric(s0))
    if (object$model$setup$include_xreg == 1) {
        k <- ncol(object$spec$xreg$xreg)
        rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
        xreg <- coredata(object$spec$xreg$xreg) %*% rho
        xreg <- xreg[start:(start + h)]
    } else {
        xreg <- rep(0, h)
    }
    xreg <- c(0, as.numeric(xreg))
    E <- matrix( rnorm(nsim*(h + 1), 0, coefficient["sigma"]), nsim, h + 1 )
    out <- simulate_aaa(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xreg)
    Y <- out$Simulated[,-1,drop = FALSE]
    if (!raw) {
        if (!is.null(object$spec$transform)) {
            Y <- object$spec$transform$inverse(as.numeric(Y), object$spec$transform$lambda)
            Y <- matrix(Y, ncol = h, nrow = nsim, byrow = FALSE)
        }
    }
    mu <- colMeans(Y)
    e <- actual - mu
    return(data.table(date = object$spec$target$index[start], error = e, horizon = 1:h))
}


hresiduals_mmm_cpp <- function(object, h = 12, nsim = 1000, start = 1, seed = NULL, raw = TRUE, ...)
{
    
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    # truncate h to the available time periods
    coefficient <- object$model$setup$parmatrix[,1]
    frequency <- object$spec$seasonal$frequency
    # last state
    N <- length(object$spec$target$y)
    if ((start + h) > N) {
        h <- N - start
    }
    model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim)
    actual <- object$spec$target$y[(start):(start + h - 1)]
    # states start one period before (because of initialization seed)
    stateinit <- object$model$states[start,,drop = FALSE]
    pars <- rep(0, 6)
    pars[1] <- stateinit[1,"Level"]
    if (model[1] == 1) pars[2] <- stateinit[1,"Trend"] else pars[2] <- 1
    pars[3] <- coefficient["alpha"]
    pars[4] <- coefficient["beta"]
    pars[5] <- coefficient["gamma"]
    pars[6] <- coefficient["phi"]
    pars <- unname(as.numeric(pars))
    if (model[2] == 1) {
        s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
    } else{ 
        s0 <- rep(1, frequency)
    }
    s0 <- unname(as.numeric(s0))
    if (object$model$setup$include_xreg == 1) {
        k <- ncol(object$spec$xreg$xreg)
        rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
        xreg <- coredata(object$spec$xreg$xreg) %*% rho
        xreg <- xreg[start:(start + h)]
    } else {
        xreg <- rep(0, h)
    }
    xreg <- c(0, as.numeric(xreg))
    E <- matrix(tsaux:::rtruncnorm(nsim*(h + 1), mu = 0, sigma = coefficient["sigma"], lb = -1), nsim, h + 1)
    out <- simulate_mmm(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xreg)
    Y <- out$Simulated[,-1,drop = FALSE]
    mu <- colMeans(Y)
    if (raw) {
        e <- actual/mu - 1
    } else{
        e <- actual - mu
    }
    return(data.table(date = object$spec$target$index[start], error = e, horizon = 1:h))
}


hresiduals_mam_cpp <- function(object, h = 12, nsim = 1000, start = 1, seed = NULL, raw = TRUE, ...)
{
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    # truncate h to the available time periods
    coefficient <- object$model$setup$parmatrix[,1]
    frequency <- object$spec$seasonal$frequency
    # last state
    N <- length(object$spec$target$y)
    if ((start + h) > N) {
        h <- N - start
    }
    model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim)
    actual <- object$spec$target$y[(start):(start + h - 1)]
    # states start one period before (because of initialization seed)
    stateinit <- object$model$states[start,,drop = FALSE]
    pars <- rep(0, 6)
    pars[1] <- stateinit[1,"Level"]
    if (model[1] == 1) pars[2] <- stateinit[1,"Trend"]
    pars[3] <- coefficient["alpha"]
    pars[4] <- coefficient["beta"]
    pars[5] <- coefficient["gamma"]
    pars[6] <- coefficient["phi"]
    pars <- unname(as.numeric(pars))
    if (model[2] == 1) {
        s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
    } else{ 
        s0 <- rep(1, frequency)
    }
    s0 <- unname(as.numeric(s0))
    if (object$model$setup$include_xreg == 1) {
        k <- ncol(object$spec$xreg$xreg)
        rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
        xreg <- coredata(object$spec$xreg$xreg) %*% rho
        xreg <- xreg[start:(start + h)]
    } else {
        xreg <- rep(0, h)
    }
    xreg <- c(0, as.numeric(xreg))
    E <- matrix(tsaux:::rtruncnorm(nsim*(h + 1), mu = 0, sigma = coefficient["sigma"], lb = -1), nsim, h + 1)
    out <- simulate_mam(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xreg)
    Y <- out$Simulated[,-1,drop = FALSE]
    mu <- colMeans(Y)
    if (raw) {
        e <- actual/mu - 1
    } else{
        e <- actual - mu
    }
    return(data.table(date = object$spec$target$index[start], error = e, horizon = 1:h))
}



hresiduals_powermam_cpp <- function(object, h = 12, nsim = 1000, start = 1, seed = NULL, raw = TRUE, ...)
{
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    # truncate h to the available time periods
    coefficient <- object$model$setup$parmatrix[,1]
    frequency <- object$spec$seasonal$frequency
    # last state
    N <- length(object$spec$target$y)
    if ((start + h) > N) {
        h <- N - start
    }
    model <- c(object$model$setup$include_trend, object$model$setup$include_seasonal, h, frequency, object$model$setup$normalized_seasonality, nsim)
    actual <- object$spec$target$y[(start):(start + h - 1)]
    # states start one period before (because of initialization seed)
    stateinit <- object$model$states[start,,drop = FALSE]
    pars <- rep(0, 8)
    pars[1] <- stateinit[1,"Level"]
    if (model[1] == 1) pars[2] <- stateinit[1,"Trend"]
    pars[3] <- coefficient["alpha"]
    pars[4] <- coefficient["beta"]
    pars[5] <- coefficient["gamma"]
    pars[6] <- coefficient["phi"]
    pars[7] <- coefficient["theta"]
    pars[8] <- coefficient["delta"]
    pars <- unname(as.numeric(pars))
    if (model[2] == 1) {
        s0 <- stateinit[1,paste0("S",0:(frequency - 1))]
    } else{ 
        s0 <- rep(1, frequency)
    }
    s0 <- unname(as.numeric(s0))
    if (object$model$setup$include_xreg == 1) {
        k <- ncol(object$spec$xreg$xreg)
        rho <- matrix(coefficient[paste0("rho",1:k)], ncol = 1)
        xreg <- coredata(object$spec$xreg$xreg) %*% rho
        xreg <- xreg[start:(start + h)]
    } else {
        xreg <- rep(0, h)
    }
    xreg <- c(0, as.numeric(xreg))
    E <- matrix(tsaux:::rtruncnorm(nsim*(h + 1), mu = 0, sigma = coefficient["sigma"], lb = -1), nsim, h + 1)
    out <- simulate_powermam(model_ = model, e_ = E, pars_ = pars, s0_ = s0, x_ = xreg)
    Y <- out$Simulated[,-1,drop = FALSE]
    mu <- colMeans(Y)
    if (raw) {
        e <- actual/mu - 1
    } else{
        e <- actual - mu
    }
    return(data.table(date = object$spec$target$index[start], error = e, horizon = 1:h))
}
