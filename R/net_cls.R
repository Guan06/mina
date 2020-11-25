################################################################################

#' Network clustering based on the sparsed adjacacency matrix.
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
#' asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
#' asv_adj <- adj(asv_norm, method = "spearman")
#' asv_cls <- net_cls(asv_adj, method = "mcl")
#' @rdname net_cls-matrix
#' @return y The cluster table.
#' @exportMethod net_cls

setMethod("net_cls", signature("matrix", "ANY", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              stop("Must specify a `method`, see `? net_cls_list`.")
          }
)

###############################################################################

#' @rdname net_cls-matrix
#' @exportMethod net_cls

setMethod("net_cls", signature("matrix", "character", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              stopifnot(
                        method %in% c("mcl", "ap"),
                        is.numeric(cutoff),
                        is.logical(neg)
              )
              if (method == "mcl" && neg == TRUE) {
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

#' Network clustering based on the sparsed adjacacency matrix.
#'
#' @include all_classes.R all_generics.R
#' @importFrom MCL mcl
#' @importFrom apcluster apcluster
#' @param x An object of class `mina` with `adj` defined.
#' @param method The clustering method used.
#' @param cutoff The cutoff for the sparsed adjacacency matrix, default 0.4.
#' @param neg Whether to keep the negative edges, cannot be TRUE when using
#' `mcl` for clustering. Default FALSE.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref", depth = 1000)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl", cutoff = 0.4, neg = FALSE)
#' maize <- net_cls(maize, method = "ap", cutoff = 0.4, neg = FALSE)
#' @return x The same `mina` class with @cls added.
#' @rdname net_cls-mina
#' @exportMethod net_cls

setMethod("net_cls", signature("mina", "ANY", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              stop("Must specify a `method`, see `? net_cls_list`.")
          }
)
###############################################################################

#' @rdname net_cls-mina
#' @exportMethod net_cls

setMethod("net_cls", signature("mina", "character", "ANY", "ANY"),
          function(x, method, cutoff = 0.4, neg = FALSE) {
              stopifnot(
                        method %in% c("mcl", "ap"),
                        is.numeric(cutoff),
                        is.logical(neg)
              )
              cls(x) <- net_cls(.adj(x), method = method, cutoff = cutoff,
                               neg = FALSE)
              return(x)
          }
)

################################################################################

#' Convert APResult (apcluster output) to dataframe.
#'
#' Modified from https://rdrr.io/github/jefferis/flycircuit/src/R/clustering.R
#' #sym-as.data.frame.APResult
#'
#' @param x an {APResult} object from \pkg{apcluster}.
#' @return y A data frame with columns `ID`, `Exemplar`, `Cluster` and
#' `Cluster_size`.
#' @keywords internal

re_format_AP <- function(x) {
    exemplars <- names(ap_exemplars(x))
    clusterids <- which(names(ap_exemplars(x)) %in% exemplars)

    clusters <- ap_clusters(x)[clusterids]
    cls <- sapply(clusters, length)
    ulc <- unlist(clusters)
    y <- data.frame(ID = names(ulc),
                    Exemplar = factor(rep(exemplars, cls)),
                    Cluster = rep(clusterids, cls),
                    Cluster_size = rep(cls, cls))
    return(y)
}

ap_exemplars <- function(x) x@exemplars
ap_clusters <- function(x) x@clusters

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
