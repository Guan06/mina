% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generics.R
\name{tina}
\alias{tina}
\title{TINA community dissimilarity used in \code{\link[mina]{com_dis}}.
Function for `tina` dissimilarity/distance calculation. Modified from Schmidt
et al., 2016.}
\usage{
tina(x, ...)
}
\arguments{
\item{x}{An matrix for dissimilarity calculation.}

\item{...}{Additional parameters.}
}
\value{
The output `tina` dissimilarity matrix.
}
\description{
TINA community dissimilarity used in \code{\link[mina]{com_dis}}.
Function for `tina` dissimilarity/distance calculation. Modified from Schmidt
et al., 2016.
}
\examples{
\dontrun{
asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
asv_dis <- com_dis(asv_norm, method = "bray")
asv_dis <- com_dis(asv_norm, method = "tina", threads = 8, nblocks = 40)
asv_tina <- tina(asv_norm, cor_method = "spearman", sim_method = "w_ja",
threads = 8, nblocks = 40)
}
}
