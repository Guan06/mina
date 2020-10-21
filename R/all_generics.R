#' @importFrom methods setClass setGeneric setMethod setRefClass as
NULL

###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, the intersect samples will be remained.
#' @param x An object of the class mina with @tab and @des defined or a
#' quantitative matrix(need parameter des in this case).
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- fit_tabs(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' }
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
#' @param depth The depth for subsampling by rarefying, 1000 by default.
#' @param replace Whether to sample with replacement (\code{TRUE} by default) or
#' without replacement (\code{FALSE}) when using method `raref`.
#' @param multi Rarefy the table for multiple times, FALSE by default, indicate
#' the times of rarefaction want to be repeated, only validate for rarefaction.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "total")
#' maize <- norm_tab(maize, method = "raref")
#' maize <- norm_tab(maize, method = "raref", depth = 1000, replace = TRUE)
#' maize <- norm_tab(maize, method = "raref", depth = 1000, replace = TRUE,
#' multi = 99)
#' }
#' @export

setGeneric("norm_tab", function(x, method, depth = 1000,
                                replace = TRUE, multi = FALSE) {
    standardGeneric("norm_tab")
})

###############################################################################

#' Calculate the adjacacency matrix of @norm by correlation.
#'
#' @param x An object of the class mina with @norm defined.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param sig (optional) The asymtotic P-values, only applicable for Pearson and
#' Spearman methods.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 500, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 200)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- adj(maize, method = "spearman", sig = FALSE)
#' maize <- adj(maize, method = "sparcc", threads = 2, nblocks = 40)
#' }
#' @export

setGeneric("adj", function(x, method, sig = FALSE, threads = 80, nblocks = 400) {
    standardGeneric("adj")
})

###############################################################################

#' Calculate the community dissimilarity / distance matrix of @norm.
#'
#' @param x An object of the class mina with @norm defined or any quantitative
#' matrix.
#' @param method The dissimilarity / distance method used, default `bray`.
#' @param threads The number of threads used for parallel running, needed for
#' method `tina`.
#' @param nblocks The number of row / column for splitted sub-matrix, needed for
#' method `tina`.
#' @examples
#' data(maize)
#' maize@tab <- maize@tab[1 : 500, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 200)
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- com_dis(maize, method = "tina", threads = 2, nblocks = 40)
#' @export

setGeneric("com_dis", function(x, method = "bray", threads = 80, nblocks = 400) {
    standardGeneric("com_dis")
})

###############################################################################

#' TINA community dissimilarity used in \code{\link[mina]{com_dis}}.
#' Function for `tina` dissimilarity/distance calculation. Modified from Schmidt
#' et al., 2016.
#'
#' @include all_classes.R all_generics.R
#' @param x An matrix for dissimilarity calculation.
#' @param cor_method The method for correlation, "pearson" and "spearman" are
#' available.
#' @param sim_method The method for similarity, "w_ja" and "uw_ja" are
#' available for weighted and unweighted Jaccard similarity respectively.
#' @param threads The number of threads used for parallel running, 80 by
#' default.
#' @param nblocks The number of row / column for splitted sub-matrix, 400 by
#' default.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 1000, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' asv_norm <- maize@norm
#' asv_tina <- tina(asv_norm, cor_method = "spearman", sim_method = "w_ja",
#' threads = 2, nblocks = 40)
#' }
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
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref", depth = 5000)
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' com_r2(maize, group = c("Compartment", "Soil", "Host_genotype"))
#' @export

setGeneric("com_r2", function(x, group) {
    standardGeneric("com_r2")
})

###############################################################################

#' Same function as `com_r2` with matrix and corresponding descriptive table as
#' input.
#' @param x Dissimilarity / distance matrix which indicate variances.
#' @param des The descriptive table of samples which define the groups.
#' @param group The name(s) of column(s) used  as experimental setup group(s) in
#' descriptive file.
#'
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' x <- maize@dis
#' des <- maize@des
#' get_r2(x, des, group = c("Compartment", "Soil"))
#' }
#' @return r2 The variance ratio cannot be explained by given groups.
#' @export

setGeneric("get_r2", function(x, des, group) {
    standardGeneric("get_r2")
})

################################################################################

#' Dimensionality reduction of community dissimilarity / distance for
#' visulization.
#'
#' @param x An object of class `mina` with @dis defined or a distance matrix.
#' @param k The dimension number after reduction.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
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
#' @param d1 The dimension be visualized in x-axis, default `1`.
#' @param d2 The dimension be visualized in y-axis, default `2`.
#' @param color The column name in @des to be used for different color groups.
#' @param shape The column name in @des to be used for different shape groups,
#' default is `NULL`.
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref", depth = 5000)
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' p1a <- com_plot(maize, match = "Sample_ID", color = "Compartment")
#' p1b <- com_plot(maize, match = "Sample_ID", d1 = 3, d2 = 4,
#' color = "Compartment")
#' p2a <- com_plot(maize, match = "Sample_ID", color = "Host_genotype")
#' p2b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 3, color =
#' "Host_genotype")
#' p3a <- com_plot(maize, match = "Sample_ID", color = "Compartment", shape =
#' "Soil")
#' p3b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 4, color =
#' "Compartment", shape = "Soil")
#' @export

setGeneric("com_plot", function(x, match, d1 = 1, d2 = 2, color, shape = NULL) {
    standardGeneric("com_plot")
})

###############################################################################

#' Visulization of components distance / dissimilarity in k dimension.
#'
#' @param x A list generated by `dmr`.
#' @param des The corresponding descriptive table.
#' @param match The column name of the components IDs in `des` with exactly the
#' same as rownames in x.
#' @param d1 The dimension be visualized in x-axis, default `1`.
#' @param d2 The dimension be visualized in y-axis, default `2`.
#' @param color The column name in `des` to be used for different color groups.
#' @param shape The column name in `des` to be used for different shape groups,
#' default `NULL`.
#' @return p The plotted figure.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' asv_dmr <- maize@dmr
#' des <- maize@des
#' p1a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment")
#' p1b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 3, d2 = 4, color =
#' "Compartment")
#' p2a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Host_genotype")
#' p2b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 3, color =
#' "Host_genotype")
#' p3a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment",
#' shape = "Soil")
#' p3b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 4, color =
#' "Compartment", shape = "Soil")
#' }
#' @export

setGeneric("pcoa_plot", function(x, des, match,
                                 d1 = 1, d2 = 2, color, shape = NULL) {
    standardGeneric("pcoa_plot")
})

################################################################################

#' Network clustering of sparsed adjacacency matrix @adj.
#'
#' @param x An object of class `mina` with @adj defined.
#' @param method The clustering method used.
#' @param cutoff The cutoff for the sparse adjacency matrix, default is 0.4.
#' @param neg Whether to keep the negative edges, default is `FALSE`.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 1000, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl")
#' maize <- net_cls(maize, method = "mcl", cutoff = 0.4, neg = FALSE)
#' maize <- net_cls(maize, method = "ap")
#' maize <- net_cls(maize, method = "ap", cutoff = 0.4, neg = FALSE)
#' }
#' @export

setGeneric("net_cls", function(x, method, cutoff = 0.4, neg = FALSE) {
    standardGeneric("net_cls")
})


################################################################################

#' Get the cluster table @cls_tab from quantitative table @norm and network
#' clustering results @cls.
#'
#' @param x_norm The normalized quantitative table used for netowrk inference
#' and clustering.
#' @param x_cls The network clustering table.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundance, default is FALSE.
#' @return x_cls The quantitative table with clusters in rows.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 1000, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' maize_norm <- maize@norm
#' maize_adj <- adj(maize_norm, method = "spearman")
#' maize_cls <- net_cls(maize_adj, method = "mcl", cutoff = 0.5)
#' maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)
#' maize_cls <- net_cls(maize_adj, method = "ap", cutoff = 0.5)
#' maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)
#' }
#' @exportMethod get_net_cls_tab

setGeneric("get_net_cls_tab", function(x_norm, x_cls, uw = FALSE) {
    standardGeneric("get_net_cls_tab")
})
################################################################################

#' Get the cluster table @cls_tab from @norm and @cls.
#'
#' @param x An object of class `mina` with @norm and @cls defined.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundances, default is FALSE.
#' @examples
#' \dontrun{
#' data(maize)
#' maize@tab <- maize@tab[1 : 500, 1 : 200]
#' maize <- norm_tab(maize, method = "raref", depth = 100)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' maize <- net_cls(maize, method = "ap", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' }
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
#' @param g_size The cutoff of group size used for filtering, default is 88.
#' @param s_size The number of samples used for network inference during
#' bootstrap and permutation (when `sig` is TRUE), it should be smaller than
#' g_size to make sure the randomness; default is 30.
#' @param rm Filtering the components present in less than `per` of the samples,
#' default is TRUE.
#' @param per The percentage of present samples for filtering, default is 0.1.
#' @param sig Whether to test the significance, skip the permutation when set as
#' FALSE, default is TRUE.
#' @param bs The times for bootstrap network inference, default is 6.
#' @param pm The times for permutation network inference, default is 6.
#' @param individual Whether to output the bootstrap and permutation results of
#' each comparison individually, default is FALSE.
#' @param out_dir The output directory if `individual` is TRUE, default is the
#' current working directory.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment")
#' }
#' @export

setGeneric("bs_pm", function(x, group, g_size = 88, s_size = 30, rm = TRUE,
                             per = 0.1, sig = TRUE, bs = 6, pm = 6,
                             individual = FALSE, out_dir = "./") {
    standardGeneric("bs_pm")
})

################################################################################

#' Calculate the network distance of @multi and test the significance when @perm
#' is defined.
#'
#' @param x An object of class `mina` with @multi (and @perm if `sig` is TRUE)
#' defined.
#' @param method The distance to be calculated, "spectra" and "Jaccard" are
#' available.
#' @param evk The first `evk` eigenvalues will be used for `spectra` distance,
#' the default is 100.
#' @param egv Wheather to output the eigenvectors for Spectral distance, the
#' defult is TRUE, only validate when `method == "spectra"`.
#' @param dir The folder to output the eigenvectors, only validate when `egv ==
#' TRUE`, default is current path.
#' @param sig Whether to test the significance, if TRUE (by default), @perm is
#' needed.
#' @param skip Whether to skip the comparison when the dimenstion of adjacency
#' matrix is smaller than setted `evk`.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment")
#' maize <- net_dis(maize, method = "spectra")
#' }
#' @export

setGeneric("net_dis", function(x, method, evk = 100, egv = TRUE, dir = "./",
                               sig = TRUE, skip = TRUE) {
    standardGeneric("net_dis")
})

################################################################################

#' Calculate the network distance of bootstrap and permutation when appliable.
#'
#' @param x The folder store the network inference results.
#' @param method The distance to be calculated, "spectra" and "Jaccard" are
#' available.
#' @param evk The first `evk` eigenvalues will be used for `spectra` distance,
#' the default is 100.
#' @param sig Whether to test the significance, if TRUE (by default),
#' permutation results should be included in the folder `x`.
#' @param skip Whether to skip the comparison when the dimenstion of adjacency
#' matrix is smaller than setted `evk`.
#' @return y The `mina` object with @dis_bs, @dis_pm and @dis_stat.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment", individual = TRUE, out_dir =
#' "./individual_bs_pm/")
#' maize_stat1 <- net_dis_indi(x = "./individual_bs_pm/", method = "spectra")
#' maize_stat2 <- net_dis_indi(x = "./individual_bs_pm/", method = "Jaccard")
#' maize_stat3 <- net_dis_indi(x = "./individual_bs_pm/", method = "spectra",
#' evk = 100, skip = TRUE)
#' }
#' @export

setGeneric("net_dis_indi", function(x, method, evk = 100,
                                    sig = TRUE, skip = TRUE) {
    standardGeneric("net_dis_indi")
})
