plot.tsets.estimate = function(x, y = NULL, ...)
{
  opar <- par()
  opar$cin <- NULL
  opar$cra <- NULL
  opar$csi <- NULL
  opar$cxy <- NULL
  opar$din <- NULL
  opar$page <- NULL
  
  d <- tsdecompose(x)
  n <- 2
  if (!is.null(d$Slope))     n <- n + 1
  if (!is.null(d$Seasonal))  n <- n + 1
  if (!is.null(d$X))         n <- n + 1

  if (x$spec$model$include_damped == 1) {
    modelx <- paste0(substr(x$spec$model$model,1,1), substr(x$spec$model$model,2,2), "d", substr(x$spec$model$model,3,3))
  } else {
    modelx <- x$spec$model$model
  }

  if (x$spec$model$type == 4) {
    modelx <- paste0("ETS [Power ", modelx,"]")
  } else{
    modelx <- paste0("ETS [", modelx,"]")
  }

  if (!is.null(x$spec$transform)) {
    lambda <-  x$spec$transform$lambda
    modelx <- paste0(modelx,"\n", "BoxCox","[lambda=",round(lambda,3),"]")
  }
  colx <- viridis_pal(option = "D", end = 0.8)(n - 1)
  layout_matrix <- matrix(1:n, nrow = n, ncol = 1)
  layout(mat = layout_matrix, heights = c(1.5, rep(1.25,n - 1)), widths = rep(2, n)) # Widths of the two columns
  par(mar = c(2,2,2,4))
  # Fitted
  plot(as.zoo(d$fitted), type = "l", ylab = "", xlab = "", col = "black", main = modelx, cex.main = 0.9, lwd = 1.5)
  lines(zoo(x$spec$target$y_orig,x$spec$target$index), col = "brown", lwd = 2, lty = 2)
  grid()
  mtext("Fitted", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  legend("topleft",c("Fitted","Actual"), col = c("black","brown"), bty = "n", lty = c(1,2), lwd = c(1,0.5), cex = 0.8)

  if (is.null(x$selection[1,"MAPE"])) {
    mape_val <- tsmetrics(x)$MAPE * 100
  } else {
    mape_val <- x$selection[1,"MAPE"] * 100
  }

  legend("bottomright",paste0("MAPE = ",round(mape_val,3),"%"), bty = "n", cex = 0.8, inset = c(0.02,.02))
  k <- 1
  # Level
  plot(as.zoo(d$Level), type = "l", ylab = "", xlab = "", col = colx[k], xaxt = "n")
  grid()
  mtext("Level", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")

  # Slope
  if (!is.null(d$Slope)) {
    k <- k + 1
    plot(as.zoo(d$Slope), type = "l", ylab = "", xlab = "", col = colx[k], xaxt = "n")
    grid()
    mtext("Slope", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  }

  # Seasonal
  if (!is.null(d$Seasonal)) {
    k <- k + 1
    plot(as.zoo(d$Seasonal), type = "l", ylab = "", xlab = "", col = colx[k], xaxt = "n")
    grid()
    mtext("Seasonal", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  }

  # External
  if (!is.null(d$X)) {
    k <- k + 1
    plot(as.zoo(d$X), type = "l", ylab = "", xlab = "", col = colx[k], xaxt = "n")
    grid()
    mtext("x-reg", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  }
  suppressWarnings(par(opar))
}

plot.tsets.simulate <- function(x, y = NULL, ...)
{
  opar <- par()
  opar$cin <- NULL
  opar$cra <- NULL
  opar$csi <- NULL
  opar$cxy <- NULL
  opar$din <- NULL
  opar$page <- NULL
  n <- 2
  if (!is.null(x$Slope))     n <- n + 1
  if (!is.null(x$Seasonal))  n <- n + 1
  colx <- viridis_pal(option = "D", end = 0.8)(n - 1)
  layout_matrix <- matrix(1:n, nrow = n, ncol = 1)
  layout(mat = layout_matrix, heights = c(1.5, rep(1.25,n - 1)), widths = rep(2, n)) # Widths of the two columns
  par(mar = c(2,2,2,4))
  # Fitted
  plot(x$Simulated)
  mtext("Simulated", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  # Level
  plot(x$Level)
  mtext("Level", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  # Slope
  if (!is.null(x$Slope)) {
    plot(x$Slope)
    mtext("Slope", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  }
  # Seasonal
  if (!is.null(x$Seasonal)) {
    plot(x$Seasonal)
    mtext("Seasonal", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  }
  suppressWarnings(par(opar))
}

plot.tsets.profile <- function(x, y = NULL, type = c("metrics", "coef"), ...)
{
  opar <- par()
  opar$cin <- NULL
  opar$cra <- NULL
  opar$csi <- NULL
  opar$cxy <- NULL
  opar$din <- NULL
  opar$page <- NULL
  type <- match.arg(type[1], c("metrics", "coef"))
  if (type == "metrics") {
    layout_matrix <- matrix(1:3, nrow = 3, ncol = 1)
    layout(mat = layout_matrix)
    par(mar = c(2,2,2,4))
    plot(x$MAPE*100, date_class = "integer", interval_quantiles = c(0.1, 0.9), main = "", xlab = "horizon")
    mtext("MAPE[%]", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    plot(x$BIAS*100, date_class = "integer", interval_quantiles = c(0.1, 0.9), main = "", xlab = "horizon")
    mtext("BIAS[%]", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    plot(x$MSLRE*100, date_class = "integer", interval_quantiles = c(0.1, 0.9), main = "", xlab = "horizon")
    mtext("MSLRE[%]", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  } else {
    Variable <- NULL
    cf <- unique(x$coef$Variable)
    cf <- cf[cf %in% c("alpha","beta","gamma","phi","theta","delta","sigma")]
    n <- length(cf)
    nf <- n2mfrow(n)
    colx <- viridis_pal(alpha = 0.5)(10)
    par(mar = c(2.5,3,2,4), mfrow = nf)
    for (i in 1:n) {
      xlim_plot <- c(min(x$coef[Variable == cf[i]]$Value, x$true.coef[cf[i]]), max(x$coef[Variable == cf[i]]$Value, x$true.coef[cf[i]]))
      xlim_dist <- (xlim_plot[2] - xlim_plot[1])/10
      xlim_plot <- c(xlim_plot[1] - xlim_dist, xlim_plot[2] + xlim_dist)
      hist(x$coef[Variable == cf[i]]$Value, main = cf[i], xlab = "", col = colx[3],  ylab = "", prob = TRUE, freq = FALSE, xaxs = "i",yaxs = "i", xlim = xlim_plot)
      abline(v = x$true.coef[cf[i]], col  = "tomato2", lty = 2, lwd = 2)
      box()
      grid()
    }
  }
  suppressWarnings(par(opar))
}