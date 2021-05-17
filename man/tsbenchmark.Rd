\name{tsbenchmark}
\alias{tsbenchmark}
\alias{tsbenchmark.tsets.spec}
\title{
Model Benchmarking
}
\description{
Benchmarks and saves the details of a model for benchmarking of speed
and unit testing.
}
\usage{
\method{tsbenchmark}{tsets.spec}(object, solver = "nlminb", 
control = list(trace = 0), autodiff = FALSE, ...)
}
\arguments{
\item{object}{
An object of class tsets.spec
}
\item{solver}{
The solver to use for estimation.
}
\item{control}{
Any additional controls passed to the solver.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsetsad package.
}
\item{\dots}{
Not currently used
}
}
\details{
This is meant mostly for internal development use.
}
\value{
A data.table
}
\author{
Alexios Galanos
}
