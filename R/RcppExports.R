# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Function for correlation coefficient calculation.
#' @param mat The input matrix for correlation calculation.
#' @return The output correlation matrix.
#' @export
cp_cor <- function(mat) {
    .Call(`_mina_cp_cor`, mat)
}

