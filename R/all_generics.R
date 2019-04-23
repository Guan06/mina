#' @param x The input mina object with quantitative tab / a matrix needed to be
#' normalized.
#' @param method The method used for the normalization of quantitative table.
#' @examples
#' norm_tab(x, method = "total")
#' @export

setGeneric("norm_tab", function(x, method, ...){
    standardGeneric("norm_tab")
})
