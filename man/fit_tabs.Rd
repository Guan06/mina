% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generics.R
\name{fit_tabs}
\alias{fit_tabs}
\title{Filter the quantitative and descriptive table to make them have the same
samples, the intersect samples will be remained.}
\usage{
fit_tabs(x)
}
\arguments{
\item{x}{An object of the class mina with `tab` and `des` defined or a
quantitative matrix(need parameter des in this case).}
}
\value{
Same `mina` object but fitted `tab` and `des` (as well as `norm` if
defined)
}
\description{
Filter the quantitative and descriptive table to make them have the same
samples, the intersect samples will be remained.
}
\examples{
data(maize)
maize <- fit_tabs(maize)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
}
