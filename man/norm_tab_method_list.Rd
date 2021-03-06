% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_tab.R
\docType{data}
\name{norm_tab_method_list}
\alias{norm_tab_method_list}
\title{List of normalization methods supported in \code{\link[mina]{norm_tab}}}
\format{
A list of character vectors.
\describe{
   \item{\code{raref}}{ By downsampling all samples to specific depth. }
   \item{\code{total}}{ Devided by the total read of each sample. }
}
}
\usage{
norm_tab_method_list
}
\description{
Normalization methods should be specified by exact string match.
}
\examples{
? norm_tab_method_list
}
\seealso{
\code{\link[mina]{norm_tab}}
}
\keyword{datasets}
