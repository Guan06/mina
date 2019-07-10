###############################################################################

#' Dimensionality reduction of the distance matrix.
#'
#' @include all_classes.R all_generics.R
#' @importFrom stats cmdscale
#' @param x A distance matrix.
#' @param k The number of dimensionality after redunction, 2 by default.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize@norm <- maize@norm[1:500, 1:300]
#' maize <- com_dis(maize, method = "bray")
#' x <- maize@dis
#' y <- dmr(x)
#' @return y The coordinates of components indicated in distance matrix in k
#' dimension.
#' @exportMethod dmr

setMethod("dmr", signature("matrix", "ANY"), function(x, k = 2) {
              y <- cmdscale(x, k = k, eig = T)
              return(y)
})

###############################################################################

#' Dimensionality reduction of the @dis included in mina.
#' @include all_classes.R all_generics.R
#'
#' @param x An object of the class `mina` with @dis defined.
#' @param k The number of dimensionality after redunction, 2 by default.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize@norm <- maize@norm[1:500, 1:300]
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' @return x The same object with @dmr added.
#' @exportMethod dmr

setMethod("dmr", signature("mina", "ANY"), function(x, k = 2) {
    x@dmr <- dmr(x@dis, k = k)
    return(x)
})
