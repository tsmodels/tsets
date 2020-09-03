\name{tsreport}
\alias{tsreport.tsets.estimate}
\alias{tsreport}
\title{
Estimation Summary pdf Report
}
\description{
Generates a pdf summary of the estimated model.
}
\usage{
\method{tsreport}{tsets.estimate}(object, output_dir = "/", args = list(name = NULL), ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsets.estimate}.
}
\item{output_dir}{
A directory where the output is written and the object saved.
}
\item{args}{
An additional list of arguments used in the generation of the report.
Only name is currently used to display the name of the data series.
}
\item{\dots}{
Not currently used.
}
}
\value{
A pdf report
}
\author{
Alexios Galanos.
}
