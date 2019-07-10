###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, the intersect samples will be remained.
#'
#' @param x An object of the class mina with @tab and @des defined or a
#' quantitative matrix(need parameter des in this case).
#' @examples
#' data(maize)
#' fit_tabs(maize)
#' @export

setGeneric("fit_tabs", function(x) {
    standardGeneric("fit_tabs")
})

###############################################################################

#' Normalize the @tab and obtain @norm for later analysis.
#'
#' @param x The input mina object with quantitative tab / a matrix needed to be
#' normalized.
#' @param method The method used for the normalization of quantitative table.
#' @examples
#' data(maize_asv)
#' maize_asv_norm <- norm_tab(maize_asv, method = "total")
#' @export

setGeneric("norm_tab", function(x, method, depth = 1000, replace = TRUE) {
    standardGeneric("norm_tab")
})

###############################################################################

#' Calculate the adjacacency matrix of @norm by correlation.
#'
#' @param x An object of the class mina with @norm defined.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- adj(maize, method = "pearson")
#' @export

setGeneric("adj", function(x, method, threads = 80, nblocks = 400) {
    standardGeneric("adj")
})

###############################################################################

#' Calculate the Person / Spearman correlation of @norm.
#'
#' @param x An matrix for pearson correlation calculation.
#' @examples
#' data(maize_asv)
#' maize_asv <- maize_asv[1:500, 1:300]
#' maize_asv_pearson <- cp_cor(maize_asv)
#' @return y The adjacacency matrix.
#' @keywords internal

setGeneric("cp_cor", function(x) {
    standardGeneric("cp_cor")
})

###############################################################################

#' Calculate the community dissimilarity / distance matrix of @norm.
#'
#' @param x An object of the class mina with @norm defined or any quantitative
#' matrix.
#' @param method The dissimilarity / distance method used.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "total")
#' maize <- com_dis(maize, method = "bray")
#' }
#' @export

setGeneric("com_dis", function(x, method, threads = 80, nblocks = 400) {
    standardGeneric("com_dis")
})

###############################################################################

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
#' \dontrun{
#' data(maize_asv)
#' maize_asv_norm <- norm_tab(maize_asv, method = "total")
#' maize_asv_tina <- tina(maize_asv_norm, cor_method = "spearman",
#' sim_method = "w_ja", threads = 80, nblocks = 400)
#' }
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
#'
#' @param x An object of class `mina` with @dis and @des defined.
#' @param group The name(s) of column(s) defined as experimental setup group(s).
#'
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "total")
#' maize <- com_dis(maize, method = "bray")
#' com_r2(maize, group = c("Compartment", "Soil", "Genotype"))
#' }
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
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "total")
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' }
#' @export

setGeneric("dmr", function(x, k = 2) {
    standardGeneric("dmr")
})

###############################################################################

#' Visulization of components distance / dissimilarity in k dimension.
#'
#' @param x An object of class `mina` with @dmr and @des defined.
#' @param match The column name of the components IDs in @des which exactly the
#' same indicated in @dmr.
#' @param color The column name in @des to be used for different color groups.
#' @param shape The column name in @des to be used for different shape groups,
#' default `NULL`.
#' shape groups.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "total")
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' p <- com_plot(maize, match = "Sample_ID", color = "Compartment")
#' p2 <- com_plot(maize, match = "Sample_ID", color = "Compartment", shape =
#' "Soil")
#' }
#' @export

setGeneric("com_plot", function(x, match, color, shape = NULL) {
    standardGeneric("com_plot")
})

################################################################################

#' Network clustering of sparsed adjacacency matrix @adj.
#'
#' @param x An object of class `mina` with @adj defined.
#' @param method The clustering method used.
#' @param cutoff The cutoff for the sparsed adjacacency matrix, default 0.4.
#' @param neg Whether to keep the negative edges, default FALSE.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl")
#' @export

setGeneric("net_cls", function(x, method, cutoff = 0.4, neg = FALSE) {
    standardGeneric("net_cls")
})

################################################################################

#' Get the cluster table @cls_tab from @tab and @cls.
#'
#' @param x An object of class `mina` with @tab and @cls defined.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundance, default FALSE.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- adj(maize, method = "spearman")
#' maize@adj <- maize@adj[1:500, 1:500]
#' maize <- net_cls(maize, method = "mcl")
#' maize@norm <- maize@norm[rownames(maize@norm) %in% maize@cls$ID, ]
#' maize <- net_cls_tab(maize)
#' @export

setGeneric("net_cls_tab", function(x, uw = FALSE) {
    standardGeneric("net_cls_tab")
})

################################################################################

#' Inferring the network of different group of samples and test significance by
#' permutation.
#'
#' @param x An object of class `mina` with @norm and @des defined.
#' @param group The column name of descriptive file @des for comparison.
#' @param g_size The cutoff of group size used for filtering, default 100.
#' @param s_size The number of samples used for network inference during
#' bootstrap and permutation (when sig == TRUE), it should be smaller than
#' g_size to make sure the randomness; default 50.
#' @param rm Filtering the components present in less than 20% of the samples,
#' default TRUE.
#' @param sig Whether to test the significance, skip the permutation when sig ==
#' FALSE, default TRUE.
#' @param bs The times for bootstrap network inference, default 6
#' @param pm The times for permuatated samples network inference, default 6.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- bs_pm(maize, group = "Compartment")
#' }
#' @export

setGeneric("bs_pm", function(x, group, g_size = 100, s_size = 50, rm = TRUE, 
                             sig = TRUE, bs = 6, pm = 6) {
    standardGeneric("bs_pm")
})

################################################################################

#' Calculate the network distance of @multi and test the significance when @perm
#' is defined.
#'
#' @param x An object of class `mina` with @multi (and @perm if sig is TRUE)
#' defined.
#' @param method The distance to be calculated, "spectral" and "jaccard" are
#' available.
#' @param sig Whether to test the significance, if TRUE (by default), @perm is
#' needed.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- bs_pm(maize, group = "Compartment")
#' maize <- net_dis(maize)
#' }
#' @export

setGeneric("net_dis", function(x, method, sig = TRUE) {
    standardGeneric("net_dis")
})
