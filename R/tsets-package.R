#' @rawNamespace useDynLib(tsets,.registration = TRUE)
#' @keywords internal
#' @import methods
#' @import tsmethods
#' @importFrom utils head tail data capture.output read.csv
#' @importFrom stats pnorm median na.omit fitted coef quantile residuals predict as.formula rpois runif rnorm acf optim qnorm sd simulate arima.sim
#' diffinv ts tsp tsp<- na.contiguous alias deviance df.residual formula hat lm model.matrix pchisq pf vcov weights
#' decompose lsfit var sigma nlminb logLik spec.ar na.pass pgamma qqline qqnorm cov
#' @importFrom graphics hist grid legend lines par plot points abline axis axis.Date axis.POSIXct box polygon layout mtext
#' @importFrom grDevices gray colorRampPalette n2mfrow
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom knitr kable
#' @importFrom data.table data.table as.data.table fwrite .N rbindlist dcast
#' @importFrom viridis viridis_pal
#' @importFrom tsaux box_cox tslinear check_xreg check_newxreg do.call.fast mape mslre mase mis bias sampling_frequency fourier_series future_dates crps tstransform auto_clean
#' @importFrom zoo index as.zoo zoo coredata na.locf
#' @importFrom xts xts as.xts is.xts endpoints
#' @importFrom stlplus stlplus
#' @importFrom corpcor make.positive.definite
#' @importFrom tsetsad estimate_ad.tsets.spec
#' @importFrom lubridate days weeks month years year %m+% tz days_in_month ymd
#' @importFrom rmarkdown render pdf_document
#' @importFrom bookdown pdf_document2
#' @importFrom Rsolnp solnp
#' @importFrom bootstrap jackknife
#' @importFrom ks kde rkde
#' @importFrom truncnorm qtruncnorm
#' @importFrom EnvStats rosnerTest
#' @importFrom Rcpp evalCpp loadModule
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
