% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj.R
\docType{data}
\name{adj_method_list}
\alias{adj_method_list}
\title{List of adjacency matix calculation methods/ orrelations supported in
\code{\link[mina]{adj}}}
\format{
A list of character vectors.
\describe{
   \item{pearson}{
       Pearson correlation.
   }
   \item{spearman}{
       Spearman correlation.
   }
   \item{sparcc}{
       SparCC correlation by spearman.
   }

}
}
\usage{
adj_method_list
}
\description{
Correlation methods should be specified by exact string match.
}
\examples{
? adj_method_list
}
\seealso{
\code{\link[mina]{adj}}
}
\keyword{datasets}
