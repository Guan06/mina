% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/com_plot.R
\name{com_plot,mina,ANY-method}
\alias{com_plot,mina,ANY-method}
\alias{com_plot,mina,character-method}
\title{Visulization of components distance / dissimilarity in k dimension.}
\usage{
\S4method{com_plot}{mina,ANY}(x, match, d1 = 1, d2 = 2, color, shape = NULL, ...)

\S4method{com_plot}{mina,character}(x, match, d1 = 1, d2 = 2, color, shape = NULL, ...)

\S4method{com_plot}{mina,character}(x, match, d1 = 1, d2 = 2, color, shape = NULL, ...)
}
\arguments{
\item{x}{An object of `mina` with list `dmr` defined.}

\item{match}{The column name of the components IDs in `des` with exactly the
same as rownames in x.}

\item{d1}{The dimension be visualized in x-axis, default `1`.}

\item{d2}{The dimension be visualized in y-axis, default `2`.}

\item{color}{The column name in `des` to be used for different color groups.}

\item{shape}{The column name in `des` to be used for different shape groups,
default `NULL`.}

\item{...}{Additional parameters.}
}
\value{
p The plotted figure.

The PCoA plot.
}
\description{
Visulization of components distance / dissimilarity in k dimension.
}
\examples{
maize <- new("mina", tab = maize_asv, des = maize_des)
maize <- norm_tab(maize, method = "raref")
maize <- fit_tabs(maize)
maize <- com_dis(maize, method = "bray")
maize <- dmr(maize)
p1a <- com_plot(maize, match = "Sample_ID", color = "Compartment")
p1b <- com_plot(maize, match = "Sample_ID", d1 = 3, d2 = 4,
color = "Compartment")
p2a <- com_plot(maize, match = "Sample_ID", color = "Host_genotype")
p2b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 3, color =
"Host_genotype")
p3a <- com_plot(maize, match = "Sample_ID", color = "Compartment", shape =
"Soil")
p3b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 4, color =
"Compartment", shape = "Soil")
}
