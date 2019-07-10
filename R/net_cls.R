################################################################################

#' Network clustering based on the sparsed adjacacency matrix @adj.
#'
#' @include all_classes.R all_generics.R
#' @importFrom MCL mcl
#' @importFrom apcluster apcluster
#' @param x Adjacency matrix used for clustering.
#' @param method The clustering method used.
#' @param cutoff The cutoff for the sparsed adjacacency matrix, default 0.4.
#' @param neg Whether to keep the negative edges, cannot be TRUE when using
#' `mcl` for clustering. Default FALSE.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize)
#' maize <- adj(maize, method = "spearman")
#' maize_adj <- maize@adj[1:500, 1:500]
#' maize_net_cls <- net_cls(maize_adj, method = "mcl")
#' @return y The cluster table.
#' @exportMethod net_cls

setMethod("net_cls", signature("matrix", "ANY", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              stop("Must specify a `method`, see `? net_cls_list`.")
          }
)

setMethod("net_cls", signature("matrix", "character", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              if (method == "mcl" && neg == "TRUE") {
                  stop("Only positive edges when using Markov Clustering.")
              }
              nr_before <- nrow(x)

              if (neg == TRUE) {
                  x[abs(x) < cutoff] <- 0
              } else {
                  x[x < cutoff] <- 0
              }

              x <- x[rowSums(x) > 0 , colSums(x) > 0]

              nr_after <- nrow(x)
              n_rm <- nr_before - nr_after
              if (n_rm > 0) {
                  message(n_rm, " components are removed for clustering.")
              }

              if (method == "mcl") {
                  mcl_cls <- mcl(x = x, addLoops = TRUE, inflation = 2.5)

                  # reformat to data frame
                  y <- re_format_MCL(mcl_cls, names = rownames(x))
              }

              if (method == "ap") {
                  # Get the degree of nodes for AP clustering.
                  x_bi <- x
                  x_bi[x_bi > 0] <- 1
                  degree <- rowSums(x_bi)
                  degree_norm <- degree / max(degree)

                  diag(x) <- 0
                  x[x == 0] <- NA

                  ap_cls <- apcluster(x, p = degree_norm)
                  y <- re_format_AP(ap_cls)
              }
              return(y)
          }
)
################################################################################

#' Network clustering based on the sparsed adjacacency matrix @adj.
#'
#' @include all_classes.R all_generics.R
#' @importFrom MCL mcl
#' @importFrom apcluster apcluster
#' @param x An object of class `mina` with @adj defined.
#' @param method The clustering method used.
#' @param cutoff The cutoff for the sparsed adjacacency matrix, default 0.4.
#' @param neg Whether to keep the negative edges, cannot be TRUE when using
#' `mcl` for clustering. Default FALSE.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- adj(maize, method = "spearman")
#' maize@adj <- maize@adj[1:500, 1:500]
#' maize <- net_cls(maize, method = "mcl")
#' @return x The same `mina` class with @cls added.
#' @exportMethod net_cls

setMethod("net_cls", signature("mina", "ANY", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              stop("Must specify a `method`, see `? net_cls_list`.")
          }
)

setMethod("net_cls", signature("mina", "character", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              x@cls <- net_cls(x@adj, method = method, cutoff = cutoff,
                               neg = FALSE)
          }
)

################################################################################

#' Convert APResult (apcluster output) to dataframe.
#'
#' Modified from https://rdrr.io/github/jefferis/flycircuit/src/R/clustering.R
#' #sym-as.data.frame.APResult
#'
#' @param x an {APResult} object from \code{\link[apcluster]{apcluster}}.
#' @return y A data frame with columns `ID`, `Exemplar`, `Cluster` and
#' `Cluster_size`.
#' @keywords internal
#' @seealso \code{\link[apcluster]{apcluster}} \code{\link[apcluster]{APResult}}

re_format_AP <- function(x) {
    exemplars <- names(x@exemplars)
    clusterids <- which(names(x@exemplars) %in% exemplars)

    clusters <- x@clusters[clusterids]
    cls <- sapply(clusters, length)
    ulc <- unlist(clusters)
    y <- data.frame(ID = names(ulc),
                    Exemplar = factor(rep(exemplars, cls)),
                    Cluster = rep(clusterids, cls),
                    Cluster_size = rep(cls, cls))
    return(y)
}

################################################################################

#' Convert mcl (mcl output) to dataframe.
#'
#' Modified from https://rdrr.io/github/jefferis/flycircuit/src/R/clustering.R
#' #sym-as.data.frame.APResult
#'
#' @param x an `mcl` object from \code{\link[MCL]{mcl}}.
#' @param names The names of clustered components.
#' @return y A data frame with columns `ID`, `Cluster` and `Cluster_size`.
#' @keywords internal
#' @seealso \code{\link[MCL]{mcl}}

re_format_MCL <- function(x, names) {
    cls <- unlist(lapply(split(x$Cluster, f = x$Cluster), length))
    y <- data.frame(ID = names,
                    Cluster = x$Cluster)
    y <- y[order(y$Cluster), ]
    y$Cluster_size = rep(cls, cls)
    return(y)
}
