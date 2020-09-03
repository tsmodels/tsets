seasonal.matrix <- function(frequency)
{
  if (frequency < 2) stop("seasonal.matrix: frequency must be > 1")
  Fmat <- matrix(0, frequency, frequency)
  Fmat[1, frequency] <- 1
  Fmat[2:frequency, 1:(frequency - 1)] <- diag(frequency - 1)
  return(Fmat)
}

tsets_filter_aaa <- function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(0, n + 1)
  smat <- matrix(0, ncol = frequency, nrow = n + 1)
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:(frequency - 1)] <- s0
    smat[1, frequency] <- -sum(s0)
  }
  err <- f <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    f[i] <- lmat[i - 1] + phi*bmat[i - 1] + smat[i - 1, frequency] + x[i]
    err[i] <- (y[i] - f[i])
    if (setup$normalized_seasonality == 1) {
      a <- gamma/frequency * err[i]
    } else{
      a <- 0
    }
    lmat[i] <- (lmat[i - 1] + phi * bmat[i - 1] + alpha * err[i]) + a
    if (setup$include_trend == 1) bmat[i] <- bmat[i - 1] * phi + beta * err[i]
    if (setup$include_seasonal == 1) {
      sx <- (smat[i - 1, frequency] + gamma * err[i]) - a
      smat[i,] <- (Fmat %*% smat[i - 1, ] ) - a
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level","Trend",paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states, residuals = err[-1]))
}


tsets_filter_mmm <- function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(1, n + 1)
  smat <- matrix(1, ncol = frequency, nrow = n + 1)
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:(frequency - 1)] <- s0
    smat[1, frequency] <- frequency - sum(s0)
  }
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  err <- f <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    f[i] <- lmat[i - 1] * bmat[i - 1] * smat[i - 1, frequency] + x[i]
    err[i] <- (y[i] - f[i])/f[i]
    if (setup$normalized_seasonality == 1) {
      a <- 1 + gamma/frequency * err[i] * smat[i - 1, frequency]
    } else {
      a <- 1
    }
    lmat[i] <- (lmat[i - 1] * (bmat[i - 1]^phi) * (1 + alpha * err[i])) * a
    if (setup$include_trend == 1) bmat[i] <- (bmat[i - 1]^phi)*(1 + beta * err[i])
    if (setup$include_seasonal == 1) {
      sx <- (smat[i - 1, frequency] * (1 + gamma * err[i])) / a
      smat[i,] <- (Fmat %*% smat[i - 1, ]) / a
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states, residuals = err[-1]))
}

tsets_filter_mam <- function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(0,n + 1)
  smat <- matrix(1, ncol = frequency, nrow = n + 1)
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:(frequency - 1)] <- s0
    smat[1, frequency] <- frequency - sum(s0)
  }
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  err <- f <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    f[i] <- (lmat[i - 1] + phi * bmat[i - 1] + x[i]) * smat[i - 1, frequency]
    err[i] <- (y[i] - f[i]) / f[i]
    if (setup$normalized_seasonality == 1) {
      a <- 1 + gamma/frequency * err[i] * smat[i - 1, frequency]
    } else {
      a <- 1
    }
    lmat[i] <- ((lmat[i - 1] + phi * bmat[i - 1]) * (1 + alpha * err[i]) ) * a
    if (setup$include_trend == 1) bmat[i] <- (phi * bmat[i - 1] + beta * (lmat[i - 1] + phi * bmat[i - 1]) * err[i]) * a
    if (setup$include_seasonal == 1) {
      sx <- (smat[i - 1, frequency] * (1 + gamma * err[i])) / a
      smat[i,] <- (Fmat %*% smat[i - 1, ]) / a
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states, residuals = err[-1]))
}

tsets_filter_powermam <- function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, theta, delta, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(0, n + 1)
  smat <- matrix(1, ncol = frequency, nrow = n + 1)
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:(frequency - 1)] <- s0
    smat[1, frequency] = frequency - sum(s0)
  }
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  err <- f <- fd <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    f[i] <- (lmat[i - 1] + phi * bmat[i - 1] + x[i]) * smat[i - 1, frequency]
    fd[i] <- ((lmat[i - 1] + phi * bmat[i - 1] + x[i])^theta) * (smat[i - 1, frequency]^delta)
    err[i] <- (y[i] - f[i]) / fd[i]
    lmat[i] <- (lmat[i - 1] + phi * bmat[i - 1]) + alpha*((lmat[i - 1] + phi * bmat[i - 1])^theta) * (smat[i - 1, frequency]^(delta - 1)) * err[i]
    if (setup$include_trend == 1) bmat[i] <- phi * bmat[i - 1] + beta * ((lmat[i - 1] + phi * bmat[i - 1])^theta) * (smat[i - 1, frequency]^(delta - 1)) * err[i]
    if (setup$include_seasonal == 1) {
      sx <- smat[i - 1, frequency] + gamma * (smat[i - 1, frequency]^delta) * ((lmat[i - 1] + phi * bmat[i - 1])^(theta - 1)) * err[i]
      smat[i,] <- Fmat %*% smat[i - 1, ]
      smat[i,1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states, residuals = err[-1], fpower = fd[-1]))
}
########################################################################################
tsfilter.tsets.estimate = function(object, y, newxreg = NULL, ...)
{
  # check what part of y is after object$spec$target$y
  yold <- xts(object$spec$target$y_orig, object$spec$target$index)
  ynew <- y
  exc <- which(index(ynew) %in% index(yold))
  if (length(exc) == 0) {
    y <- ynew
  } else {
    y <- ynew[-exc]
    if (NROW(y) == 0) {
      warning("\nno new data in y which is not already in object!")
      return(object)
    }
  }
  newindex <- index(y)
  yneworig <- y
  if (!is.null(object$spec$transform)) {
    y <- object$spec$transform$transform(yneworig, object$spec$transform$lambda)
  } else{
    y <- yneworig
  }
  alpha <- object$model$setup$parmatrix["alpha",1]
  beta <- object$model$setup$parmatrix["beta",1]
  gamma <- object$model$setup$parmatrix["gamma",1]
  phi <- object$model$setup$parmatrix["phi",1]
  delta <- object$model$setup$parmatrix["delta",1]
  theta <- object$model$setup$parmatrix["theta",1]
  frequency <- object$spec$seasonal$frequency
  setup <- object$model$setup
  pnames <- rownames(object$model$setup$parmatrix)
  s0 <- as.numeric(tail(object$model$states[,paste0("S",0:(frequency - 1))],1))
  l0 <- tail(object$model$states[,"Level"],1)
  b0 <- tail(object$model$states[,"Trend"],1)
  if (!is.null(newxreg) & setup$include_xreg == 1) {
    newxreg <- newxreg[index(y)]
    if (NROW(y) != NROW(newxreg)) stop("\ny and newxreg do not have the same number of rows")
  }
  if (setup$include_xreg == 1) {
    if (!is.null(newxreg)) {
      x <- newxreg %*% object$model$setup$parmatrix[grepl("[rho][0-9]",pnames),1,drop = FALSE]
    } else {
      x <- rep(0, NROW(y))
    }
  } else {
    x <- rep(0, NROW(y))
  }
  x <- c(0, x)
  y <- as.numeric(y)
  x <- as.numeric(x)
  out <- switch(object$model$setup$type,
               "1" = onlinefilter_aaa(y, alpha = alpha, beta = beta, gamma = gamma, phi = phi, l0 = l0, b0 = b0, s0 = s0, frequency = frequency, x = x, setup = setup),
               "2" = onlinefilter_mmm(y, alpha = alpha, beta = beta, gamma = gamma, phi = phi, l0 = l0, b0 = b0, s0 = s0, frequency = frequency, x = x, setup = setup),
               "3" = onlinefilter_mam(y, alpha = alpha, beta = beta, gamma = gamma, phi = phi, l0 = l0, b0 = b0, s0 = s0, frequency = frequency, x = x, setup = setup),
               "4" = onlinefilter_powermam(y, alpha = alpha, beta = beta, gamma = gamma, phi = phi, l0 = l0, b0 = b0, s0 = s0, frequency = frequency,
                                           theta = theta, delta = delta, x = x, setup = setup))
  # update the whole object to incorporate the new information
  object$spec$target$y_orig <- c(object$spec$target$y_orig, as.numeric(yneworig))
  object$spec$target$index <- c(object$spec$target$index, newindex)
  object$spec$target$y <- c(object$spec$target$y, as.numeric(y))
  if (!is.null(newxreg)) {
    object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, coredata(newxreg))
    object$spec$xreg$index <- c(object$spec$xreg$index, index(newxreg))
  }
  object$model$filtered <- c(object$model$filtered,out$filtered)
  object$model$states <- rbind(object$model$states, out$states)
  object$model$residuals <- c(object$model$residuals,out$residuals)
  return(object)
}

onlinefilter_aaa = function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(0, n + 1)
  smat <- matrix(0, ncol = frequency, nrow = n + 1)
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:frequency] <- s0
  }
  err <- f <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    if (setup$normalized_seasonality) {
      a <- gamma/frequency * err[i]
    } else {
      a <- 0
    }
    f[i] <- lmat[i - 1] + phi * bmat[i - 1] + smat[i - 1, frequency] + x[i]
    err[i] <- (y[i] - f[i])
    lmat[i] <- (lmat[i - 1] + phi * bmat[i - 1] + alpha * err[i]) + a
    if (setup$include_trend == 1) bmat[i] <- bmat[i - 1] * phi + beta * err[i]
    if (setup$include_seasonal == 1) {
      sx <- (smat[i - 1, frequency] + gamma * err[i]) - a
      smat[i,] <- (Fmat %*% smat[i - 1,]) - a
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states[-1,,drop = FALSE], residuals = err[-1]))
}


onlinefilter_mmm = function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(1, n + 1)
  smat <- matrix(1, ncol = frequency, nrow = n + 1)
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:frequency] <- s0
  }
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  err <- f <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    f[i] <- lmat[i - 1] * (bmat[i - 1]^phi) * smat[i - 1, frequency] + x[i]
    err[i] <- (y[i] - f[i]) / f[i]
    if (setup$normalized_seasonality) {
      a <- 1 + gamma/frequency * err[i] * smat[i - 1, frequency]
    } else {
      a <- 1
    }
    lmat[i] <- (lmat[i - 1] * (bmat[i - 1]^phi) * (1 + alpha * err[i])) * a
    if (setup$include_trend == 1) bmat[i] <- (bmat[i - 1]^phi) * (1 + beta * err[i])
    if (setup$include_seasonal == 1) {
      sx <- (smat[i - 1, frequency] * (1 + gamma * err[i]) ) / a
      smat[i,] <- (Fmat %*% smat[i - 1, ]) / a
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states[-1, , drop = FALSE], residuals = err[-1]))
}

onlinefilter_mam = function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(0, n + 1)
  smat <- matrix(1, ncol = frequency, nrow = n + 1)
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:frequency] <- s0
  }
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0

  err <- f <- rep(0, n + 1)

  y <- c(0,as.numeric(y))

  for (i in 2:(n + 1)) {
    f[i] <- (lmat[i - 1] + phi * bmat[i - 1] + x[i]) * smat[i - 1, frequency]
    err[i] <- (y[i] - f[i]) / f[i]
    if (setup$normalized_seasonality) {
      a <- 1 + gamma/frequency * err[i] * smat[i - 1, frequency]
    } else {
      a <- 1
    }
    lmat[i] <- ((lmat[i - 1] + phi * bmat[i - 1]) * (1 + alpha * err[i])) * a
    if (setup$include_trend == 1) bmat[i] <- (phi * bmat[i - 1] + beta * (lmat[i - 1] + phi * bmat[i - 1]) * err[i]) * a
    if (setup$include_seasonal == 1) {
      sx <- (smat[i - 1, frequency] * (1 + gamma * err[i])) / a
      smat[i,] <- (Fmat %*% smat[i - 1,]) / a
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states[-1, , drop = FALSE], residuals = err[-1]))
}

onlinefilter_powermam = function(y, alpha, beta, gamma, phi = 1, l0, b0, s0, frequency = 12, theta, delta, x, setup)
{
  n <- length(y)
  lmat <- bmat <- rep(0, n + 1)
  smat <- matrix(1, ncol = frequency, nrow = n + 1)
  if (setup$include_seasonal == 1) {
    Fmat <- seasonal.matrix(frequency)
    smat[1, 1:frequency] <- s0
  }
  lmat[1] <- l0
  if (setup$include_trend == 1) bmat[1] <- b0
  err <- f <- fd <- rep(0, n + 1)
  y <- c(0, as.numeric(y))
  for (i in 2:(n + 1)) {
    f[i] <- (lmat[i - 1] + phi * bmat[i - 1] + x[i]) * smat[i - 1, frequency]
    fd[i] <- ((lmat[i - 1] + phi * bmat[i - 1] + x[i])^theta) * (smat[i - 1, frequency]^delta)
    err[i] <- (y[i] - f[i]) / fd[i]
    lmat[i] <- (lmat[i - 1] + phi * bmat[i - 1]) + alpha * ((lmat[i - 1] + phi * bmat[i - 1])^theta) * (smat[i - 1, frequency]^(delta - 1)) * err[i]
    if (setup$include_trend == 1) bmat[i] = phi * bmat[i - 1] + beta * ((lmat[i - 1] + phi * bmat[i - 1])^theta) * (smat[i - 1, frequency]^(delta - 1)) * err[i]
    if (setup$include_seasonal == 1) {
      sx <- smat[i - 1, frequency] + gamma * (smat[i - 1, frequency]^delta) * ((lmat[i - 1] + phi * bmat[i - 1])^(theta - 1)) * err[i]
      smat[i,] <- Fmat %*% smat[i - 1,]
      smat[i, 1] <- sx
    }
  }
  states <- cbind(lmat, bmat, smat)
  colnames(states) <- c("Level", "Trend", paste0("S",0:(frequency - 1)))
  return(list(filtered = f[-1], states = states[-1, , drop = FALSE], residuals = err[-1], fpower = fd[-1]))
}
