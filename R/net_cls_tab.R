################################################################################

#' Get the cluster table @cls_tab from quantitative table @tab and network
#' clustering results @cls.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class `mina` with @norm_tab and @cls defined.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundance, default FALSE.
#' @return x The same `mina` object with @cls_tab added.
#' @examples
#' data(maize)
#' maize@tab <- maize@tab[1 : 1000, 1 : 500]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' maize <- net_cls(maize, method = "ap", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' @exportMethod net_cls_tab

setMethod("net_cls_tab", signature("mina", "ANY"),
          function(x, uw = FALSE) {
              tab <- x@norm
              cls <- x@cls

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
              x@cls_tab <- cls_tab
              return(x)
          }
)
