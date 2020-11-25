###############################################################################

#' Dimensionality reduction of the distance matrix.
#'
#' @include all_classes.R all_generics.R
#' @importFrom stats cmdscale
#' @param x A distance matrix.
#' @param k The number of dimensionality after reduction, 4 by default.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' asv_dis <- dis(maize)
#' asv_dis_dmr <- dmr(asv_dis, k = 4)
#' @return y The coordinates of components indicated in distance matrix in k
#' dimension.
#' @rdname dmr-matrix
#' @exportMethod dmr

setMethod("dmr", signature("matrix", "ANY"), function(x, k = 4) {
    stopifnot(is.numeric(k))
    x[is.na(x)] <- 0
    y <- cmdscale(x, k = k, eig = T)
    return(y)
})

###############################################################################

#' Dimensionality reduction of the `dis` included in mina.
#' @include all_classes.R all_generics.R
#'
#' @param x An object of the class `mina` with `dis` defined.
#' @param k The number of dimensionality after reduction, 4 by default.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' @return x The same object with `dmr` added.
#' @rdname dmr-mina
#' @exportMethod dmr

setMethod("dmr", signature("mina", "ANY"), function(x, k = 4) {
    stopifnot(is.numeric(k))
    .dmr(x) <- dmr(dis(x), k = k)
    return(x)
})
