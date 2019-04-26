###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, the intersect samples will be remained.
#'
#' @param x An object of the class mina with @tab and @des_tab defined or a
#' quantitative matrix(need parameter des in this case).
#' @examples
#' fit_tabs(x)
#' @export

setGeneric("fit_tabs", function(x) {
    standardGeneric("fit_tabs")
})

###############################################################################

#' @param x The input mina object with quantitative tab / a matrix needed to be
#' normalized.
#'
#' @param method The method used for the normalization of quantitative table.
#' @examples
#' norm_tab(x, method = "total")
#' @export

setGeneric("norm_tab", function(x, method, depth = 1000, replace = TRUE) {
    standardGeneric("norm_tab")
})
