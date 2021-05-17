tsbenchmark.tsets.spec <- function(object, solver = "nlminb", control = list(trace = 0), autodiff = FALSE, ...)
{
    start <- Sys.time()
    mod <- estimate(object, solver = solver, control = control, autodiff = autodiff)
    end <- Sys.time()
    return(data.table(start = start, end = end, spec = list(object), estimate = list(mod), solver = solver, control = control, loglik = as.numeric(logLik(mod))))
}
