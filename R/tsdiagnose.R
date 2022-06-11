check_parameters <- function(object)
{
    error_type <- object$spec$model$error
    # Additive
    cf <- coef(object)
    cf_names <- names(cf)
    r <- residuals(object, raw = TRUE)
    r <- as.numeric(na.omit(r))
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
    ipars <- object$model$setup$parmatrix
    lb_check <- ipars[, "init"] >= ipars[, "lower"]
    lb_check <- lb_check[c("alpha","beta","gamma","phi","theta","delta")]
    ub_check <- ipars[, "init"] <= ipars[, "upper"]
    ub_check <- ub_check[c("alpha","beta","gamma","phi","theta","delta")]

    if (error_type == "Additive") {
        condition_table <- data.table(coef = c("alpha","beta","gamma","phi","theta","delta"), value = round(as.numeric(cf[c("alpha","beta","gamma","phi","theta","delta")]),4),
                   ">lb" = lb_check, "<ub" = ub_check, "condition" = c("NA"," < alpha"," < (1 - alpha)","NA","NA","NA"), "condition_pass" = c(NA, condition_slope, condition_seasonal,NA, NA, NA))
    } else {
        condition_residuals <- all(r > (-1))
        condition_table <- data.table(coef = c("alpha","beta","gamma","phi","theta","delta","residuals"), value = c(round(as.numeric(cf[c("alpha","beta","gamma","phi","theta","delta")]),4), NA),
                   ">lb" = c(lb_check, condition_residuals), "<ub" = c(ub_check,NA))
    }
    return(condition_table)
}

#' Model Diagnostics
#'
#' @description Creates a short summary of model based diagnostics.
#' @param object an object of class \dQuote{tsets.estimate}.
#' @param plot whether to generate diagnostic plots to accompany summary.
#' @param ... not currently used.
#' @return A list of tables (printed out and returned invisibly) with
#' Ljung-Box test for residual autocorrelation, parameter and model bounds
#' diagnostics and outlier dates using the Rosner test from the
#' \code{EnvStats} package.
#' @aliases tsdiagnose
#' @method tsdiagnose tsets.estimate
#' @rdname tsdiagnose
#' @export
#'
tsdiagnose.tsets.estimate <- function(object, plot = FALSE, ...)
{
    ctable <- check_parameters(object)
    r <- scale(as.numeric(na.omit(residuals(object))), center = FALSE, scale = TRUE)
    j <- 0
    df <- 0
    b1 <- weighted_box_test(r, lag = 1, type = "Ljung-Box", fitdf = 0)
    b2 <- weighted_box_test(r, lag = 2 + df + j, type = "Ljung-Box", fitdf = df)
    b3 <- weighted_box_test(r, lag = 3 + df + j, type = "Ljung-Box", fitdf = df)
    b4 <- weighted_box_test(r, lag = 4 + df + j, type = "Ljung-Box", fitdf = df)
    lbsr <- data.table(Lag =  c("[1]", "[2]", "[3]", "[4]"), statistic = c(b1$statistic[[1]], b2$statistic[[1]], b3$statistic[[1]],b4$statistic[[1]]),
                       pvalue = c(b1$p.value[[1]], b2$p.value[[1]], b3$p.value[[1]], b4$p.value[[1]]))
    cat("\nWeighted Ljung-Box Test [scaled residuals]")
    cat("\n------------------------------------------\n")
    print(lbsr, row.names = FALSE, digits = 3)
    cat("\nParameter Bounds and Conditions")
    cat("\n------------------------------------------\n")
    print(ctable, row.names = FALSE, digits = 3)
    if (object$spec$model$error != "Additive") {
        cat("\n(raw model residuals lower bound = -1)\n")
    }
    rtest <- rosnerTest(as.numeric(na.omit(residuals(object))), k = 10)
    if (any(rtest$all.stats$Outlier)) {
        out.index <- object$spec$target$index[which(object$model$setup$good == 1)][rtest$all.stats$Obs.Num[rtest$all.stats$Outlier]]
        cat("\nOutlier Diagnostics (based on Rosner Test)")
        cat("\n------------------------------------------")
        cat("\nOutliers:", as.character(out.index))
    } else {
        out.index <- NULL
    }
    if (plot) {
        rw <- na.omit(residuals(object, raw = TRUE))
        par(mfrow = c(3,1), mar = c(3,3,3,3))
        acf(as.numeric(r), type = "correlation", main = "Residuals Autocorrelation")
        hist(rw, breaks = "fd", main = "Model Residuals Histogram", probability = T)
        box()
        qqnorm(rw, main = "Model Residuals Normal Q-Q Plot")
        qqline(rw, col = 2)
    }
    L <- list(lb_test = lbsr, outliers = rtest$all.stats, outlier_index = out.index, parameter_bounds = ctable)
    return(invisible(L))
}

##########################################################################################
# Direct Import of Weighted Tests of FISHER and GALLAGHER (WeightedPortTest package)
weighted_box_test = function(x, lag = 1, type = c("Box-Pierce", "Ljung-Box", "Monti"),
                             fitdf = 0, sqrd.res = FALSE, log.sqrd.res = FALSE, abs.res = FALSE,
                             weighted = TRUE)
{
    if (lag < (2 * fitdf + fitdf - 1)) stop("\nLag must be equal to a minimum of 2*fitdf+fitdf-1")
    if (NCOL(x) > 1) stop("\nx is not a vector or univariate time series")
    if (lag < 1) stop("\nLag must be positive")
    if (fitdf < 0) stop("\nFitdf cannot be negative")
    if ((sqrd.res && log.sqrd.res) || (sqrd.res && abs.res) || (log.sqrd.res && abs.res)) stop("Only one option of: sqrd.res, log.sqrd.res or abs.res can be selected")
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)
    if (abs.res) {
        x <- abs(x)
    }
    if (sqrd.res || log.sqrd.res) {
        x <- x^2
    }
    if (log.sqrd.res) {
        x <- log(x)
    }
    if (weighted) {
        if (type == "Monti") {
            METHOD <- "Weighted Monti test (Gamma Approximation)"
            cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE,
                       na.action = na.pass)
            obs <- cor$acf[1:lag]
        }
        else {
            cor <- acf(x, lag.max = lag, type = "correlation",
                       plot = FALSE, na.action = na.pass)
            obs <- cor$acf[2:(lag + 1)]
        }
        if (type == "Ljung-Box") {
            METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
        }
        n <- sum(!is.na(x))
        weights <- (lag - 1:lag + 1)/(lag)
        if (type == "Box-Pierce") {
            METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
            STATISTIC <- n * sum(weights * obs^2)
        }
        else {
            STATISTIC <- n * (n + 2) * sum(weights * (1/seq.int(n - 1, n - lag) * obs^2))
        }
        if (sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
        }
        else if (log.sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "Weighted X-squared on Log-Squared Residuals for detecting nonlinear processes"
        }
        else if (abs.res) {
            fitdf <- 0
            names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes"
        }
        else {
            names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
        }
        shape <- (3/4) * (lag + 1)^2 * lag/(2 * lag^2 + 3 * lag + 1 - 6 * lag * fitdf)
        scale <- (2/3) * (2 * lag^2 + 3*lag + 1 - 6 * lag * fitdf)/(lag*(lag + 1))
        PARAMETER <- c(shape, scale)
        names(PARAMETER) <- c("Shape", "Scale")
        PVAL <- 1 - pgamma(STATISTIC, shape = shape, scale = scale)
        names(PVAL) <- "Approximate p-value"
    }
    else {
        if (type == "Monti") {
            METHOD <- "Monti test"
            cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE,
                       na.action = na.pass)
            obs <- cor$acf[1:lag]
        }
        else {
            cor <- acf(x, lag.max = lag, type = "correlation",
                       plot = FALSE, na.action = na.pass)
            obs <- cor$acf[2:(lag + 1)]
        }
        if (type == "Ljung-Box") {
            METHOD <- "Ljung-Box test"
        }
        n <- sum(!is.na(x))
        if (type == "Box-Pierce") {
            METHOD <- "Box-Pierce test"
            STATISTIC <- n * sum(obs^2)
        }
        else {
            STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1,
                                                      n - lag) * obs^2))
        }
        if (sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
        }
        else if (log.sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "X-squared on Log-Squared Residuals for detecting nonlinear processes"
        }
        else if (abs.res) {
            fitdf <- 0
            names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes"
        }
        else {
            names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
        }
        mydf <- lag - fitdf
        PARAMETER <- c(mydf)
        names(PARAMETER) <- c("df")
        PVAL <- 1 - pchisq(STATISTIC, df = mydf)
        names(PVAL) <- "p-value"
    }
    structure(list(statistic = STATISTIC, parameter = PARAMETER,
                   p.value = PVAL, method = METHOD, data.name = DNAME),
              class = "htest")
}
