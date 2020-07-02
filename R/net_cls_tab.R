################################################################################

#' Get the cluster table @cls_tab from quantitative table @norm and network
#' clustering results @cls.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class `mina` with @norm and @cls defined.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundance, default is FALSE.
#' @return x The same `mina` object with @cls_tab added.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 1000, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' maize <- net_cls(maize, method = "ap", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' }
#' @exportMethod net_cls_tab

setMethod("net_cls_tab", signature("mina", "ANY"),
          function(x, uw = FALSE) {
              x@cls_tab <- get_net_cls_tab(x@norm, x@cls, uw)
              return(x)
          }
)
################################################################################

#' Get the cluster table @cls_tab from quantitative table @norm and network
#' clustering results @cls.
#'
#' @include all_classes.R all_generics.R
#' @param x_norm The normalized quantitative table used for netowrk inference
#' and clustering.
#' @param x_cls The network clustering table.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundance, default is FALSE.
#' @return x_cls The quantitative table with clusters in rows.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 1000, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' maize_norm <- maize@norm
#' maize_adj <- adj(maize_norm, method = "spearman")
#' maize_cls <- net_cls(maize_adj, method = "mcl", cutoff = 0.5)
#' maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)
#' maize_cls <- net_cls(maize_adj, method = "ap", cutoff = 0.5)
#' maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)
#' }
#' @exportMethod get_net_cls_tab

setMethod("get_net_cls_tab", signature("matrix", "data.frame", "ANY"),
          function(x_norm, x_cls, uw = FALSE) {
              tab <- x_norm
              cls <- x_cls

              if (nrow(tab) != nrow(cls)) {
                  message("Different number of components!")
                  message("Filtering...")
                  tab <- tab[rownames(tab) %in% cls$ID, ]
              }

              if (uw) {
                  tab[tab > 0] <- 1
              }
              group <- cls$Cluster
              tab <- tab[match(cls$ID, rownames(tab)), ]
              cls_tab <- apply(tab, 2, function(x) rowsum(x, group))
              rownames(cls_tab) <- paste0("Cluster_", unique(group))
              return(cls_tab)
          }
)
