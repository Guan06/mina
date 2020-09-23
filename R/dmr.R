###############################################################################

#' Dimensionality reduction of the distance matrix.
#'
#' @include all_classes.R all_generics.R
#' @importFrom stats cmdscale
#' @param x A distance matrix.
#' @param k The number of dimensionality after reduction, 4 by default.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' asv_dis <- maize@dis
#' asv_dis_dmr <- dmr(asv_dis)
#' asv_dis_dmr <- dmr(asv_dis, k = 4)
#' }
#' @return y The coordinates of components indicated in distance matrix in k
#' dimension.
#' @rdname dmr-matrix
#' @exportMethod dmr

setMethod("dmr", signature("matrix", "ANY"), function(x, k = 4) {
              x[is.na(x)] <- 0
              y <- cmdscale(x, k = k, eig = T)
              return(y)
})

###############################################################################

#' Dimensionality reduction of the @dis included in mina.
#' @include all_classes.R all_generics.R
#'
#' @param x An object of the class `mina` with @dis defined.
#' @param k The number of dimensionality after reduction, 4 by default.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' maize <- dmr(maize, k = 4)
#' }
#' @return x The same object with @dmr added.
#' @rdname dmr-mina
#' @exportMethod dmr

setMethod("dmr", signature("mina", "ANY"), function(x, k = 4) {
    x@dmr <- dmr(x@dis, k = k)
    return(x)
})
