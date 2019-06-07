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

#' TINA calculation used in \code{\link[mina]{com_dis}}.
#'
#' @param x An matrix for `tina` dissimilarity calculation.
#' @param cor_method The method for correlation, "pearson" and "spearman" are
#' available.
#' @param sim_method The method for similarity, "w_ja" and "uw_ja" are
#' available for weighted and unweighted Jaccard similarity respectively.
#' @param threads (optional) The number of threads used for parallel running,
#' 80 by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' t <- tina(x, cor_method = "spearman", sim_method = "w_ja", threads = 80,
#'           nblocks = 400)
#' @return t The output `tina` dissimilarity matrix.
#' @export
setGeneric("tina", function(x, cor_method = "spearman", sim_method = "w_ja",
                            threads = 80, nblocks = 400) {
    standardGeneric("tina")
})

################################################################################

#' Calculate the unexplained variance ratio using formula indicated in:
#' Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of
#' variance. Austral Ecology, 26: 32--46.

#' @param x An object of class `mina` with @dis and @des_tab defined.
#' @param group The name(s) of column(s) defined as experimental setup group(s).
#'
#' @examples
#' x <- com_r2(x, group = "soil")
#' @export

setGeneric("com_r2", function(x, group) {
    standardGeneric("com_r2")
})

################################################################################

#' Dimensionality reduction of community dissimilarity / distance for
#' visulization.
#'
#' @param x An object of class `mina` with @dis defined.
#' @param k The dimension number after reduction.
#' @examples
#' x <- dmr(x, k = 2)
#' @export

setGeneric("dmr", function(x, k = 2) {
    standardGeneric("dmr")
})
