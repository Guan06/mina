#' @param mina A mina object including Tab.
#' @param method The method used for the normalization of quantitative table.
#' @examples
#' mina <- normTab(mina, method="raref")
#' @importFrom phyloseq rarefy_even_depth phyloseq_transform_css
#' @export

setGeneric("normTab", function(mina, method="raref", ...) standardGeneric("normTab"))
