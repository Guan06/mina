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

#' Normalize the @tab and obtain @norm for later analysis.
#' @param x The input mina object with quantitative tab / a matrix needed to be
#' normalized.
#' @param method The method used for the normalization of quantitative table.
#' @examples
#' x <- norm_tab(x, method = "total")
#' @export

setGeneric("norm_tab", function(x, method, depth = 1000, replace = TRUE) {
    standardGeneric("norm_tab")
})

###############################################################################

#' Calculate the adjacacency matrix of @norm by correlation.
#' @param x An object of the class mina with @norm defined or any quantitative
#' matrix.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' x <- adj(x, method = "pearson")
#' @export

setGeneric("adj", function(x, method, threads = 80, nblocks = 400) {
    standardGeneric("adj")
})

###############################################################################

#' Calculate the community dissimilarity / distance matrix of @norm.
#' @param x An object of the class mina with @norm defined or any quantitative
#' matrix.
#' @param method The dissimilarity / distance method used.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' x <- adj(x, method = "bray")
#' @export

setGeneric("com_dis", function(x, method, threads = 80, nblocks = 400) {
    standardGeneric("com_dis")
})
