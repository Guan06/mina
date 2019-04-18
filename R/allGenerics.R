#' @param mina A mina object including Tab.
#' @param method The method used for the normalization of quantitative table.
#' @examples
#' norm_tab(mina, "total")
#' @export

setGeneric("norm_tab", function(mina, method, ...){
    standardGeneric("norm_tab")
})
