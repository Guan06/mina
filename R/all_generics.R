#' @importFrom methods setClass setGeneric setMethod setRefClass as
NULL

###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, the intersect samples will be remained.
#' @param x An object of the class mina with `tab` and `des` defined or a
#' quantitative matrix(need parameter des in this case).
#' @return Same `mina` object but fitted `tab` and `des` (as well as `norm` if
#' defined)
#' @examples
#' data(maize)
#' maize <- fit_tabs(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' @export

setGeneric("fit_tabs", function(x) {
    standardGeneric("fit_tabs")
})

###############################################################################

#' Normalize the slot `tab` for later analysis.
#'
#' @param x The input mina object with quantitative tab / a matrix needed to be
#' normalized.
#' @param method The method used for the normalization of quantitative table.
#' @param ... Additional parameters.
#' @return Normalized quantitative table.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "total")
#' @export

setGeneric("norm_tab", function(x, method, ...) {
    standardGeneric("norm_tab")
})

###############################################################################

#' Calculate the correlation adjacacency matrix.
#'
#' @param x An object of the class mina with `norm` defined or a `norm` matrix.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param ... Additional parameters.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref", depth = 1000)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman", sig = FALSE)
#' @return Adjacency matrix.
#' @export

setGeneric("adj", function(x, method, ...) {
    standardGeneric("adj")
})

###############################################################################

#' Calculate the community dissimilarity / distance matrix.
#'
#' @param x An object of the class mina with `norm` defined or any quantitative
#' matrix.
#' @param method The dissimilarity / distance method used, default `bray`.
#' @param ... Additional parameters.
#' @examples
#' asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
#' asv_dis <- com_dis(asv_norm, method = "bray")
#' @return The distance / dissimilarity matrix.
#' @export

setGeneric("com_dis", function(x, method = "bray", ...) {
    standardGeneric("com_dis")
})

###############################################################################

#' TINA community dissimilarity used in \code{\link[mina]{com_dis}}.
#' Function for `tina` dissimilarity/distance calculation. Modified from Schmidt
#' et al., 2016.
#'
#' @include all_classes.R all_generics.R
#' @param x An matrix for dissimilarity calculation.
#' @param ... Additional parameters.
#' @examples
#' \dontrun{
#' asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
#' asv_dis <- com_dis(asv_norm, method = "bray")
#' asv_dis <- com_dis(asv_norm, method = "tina", threads = 8, nblocks = 40)
#' asv_tina <- tina(asv_norm, cor_method = "spearman", sim_method = "w_ja",
#' threads = 8, nblocks = 40)
#' }
#' @return The output `tina` dissimilarity matrix.
#' @export

setGeneric("tina", function(x, ...) {
    standardGeneric("tina")
})

################################################################################

#' Calculate the unexplained variance ratio using formula indicated in:
#' Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of
#' variance. Austral Ecology, 26: 32--46.
#'
#' @param x An object of class `mina` with `dis` and `des` defined.
#' @param group The name(s) of column(s) defined as experimental setup group(s).
#' @examples
#' data(maize)
#' maize <- norm_tab(maize, method = "raref", depth = 5000)
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' com_r2(maize, group = c("Compartment", "Soil", "Host_genotype"))
#' @return Unexplained variance ratio.
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
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' get_r2(dis(maize), des(maize), group = c("Compartment", "Soil"))
#' @return r2 The variance ratio cannot be explained by given groups.
#' @export

setGeneric("get_r2", function(x, des, group) {
    standardGeneric("get_r2")
})

################################################################################

#' Dimensionality reduction of community dissimilarity / distance for
#' visulization.
#'
#' @param x An object of class `mina` with `dis` defined or a distance matrix.
#' @param k The dimension number after reduction.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' @return The dimentionality reduction results.
#' @export

setGeneric("dmr", function(x, k = 2) {
    standardGeneric("dmr")
})

###############################################################################

#' Visulization of components distance / dissimilarity in k dimension.
#'
#' @param x An object of class `mina` with `dmr` and `des` defined.
#' @param match The column name of the components IDs in `des` which exactly the
#' same indicated in `dmr`.
#' @param ... Additional parameters.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
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
#' @return The PCoA plot.
#' @export

setGeneric("com_plot", function(x, match, ...) {
    standardGeneric("com_plot")
})

###############################################################################

#' Visulization of components distance / dissimilarity in k dimension.
#'
#' @param x A list generated by `dmr`.
#' @param des The corresponding descriptive table.
#' @param match The column name of the components IDs in `des` with exactly the
#' same as rownames in x.
#' @param ... Additional parameters.
#' @return p The plotted figure.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' asv_dmr <- .dmr(maize)
#' des <- des(maize)
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
#' @export

setGeneric("pcoa_plot", function(x, des, match, ...) {
    standardGeneric("pcoa_plot")
})

################################################################################

#' Network clustering of sparsed adjacacency matrix.
#'
#' @param x An object of class `mina` with `adj` defined.
#' @param method The clustering method used.
#' @param ... Additional parameters.
#' @return The network clustering results.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref", depth = 1000)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "mcl", cutoff = 0.4, neg = FALSE)
#' @export

setGeneric("net_cls", function(x, method, ...) {
    standardGeneric("net_cls")
})

################################################################################

#' Get the cluster table `cls_tab` from quantitative table `norm` and network
#' clustering results `cls`.
#'
#' @param x_norm The normalized quantitative table used for netowrk inference
#' and clustering.
#' @param x_cls The network clustering table.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundance, default is FALSE.
#' @return x_cls The quantitative table with clusters in rows.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref", depth = 1000)
#' maize <- fit_tabs(maize)
#' maize_norm <- norm(maize)
#' maize_adj <- adj(maize_norm, method = "spearman")
#' maize_cls <- net_cls(maize_adj, method = "ap", cutoff = 0.5)
#' maize_cls_tab <- get_net_cls_tab(maize_norm, maize_cls)
#' @exportMethod get_net_cls_tab

setGeneric("get_net_cls_tab", function(x_norm, x_cls, uw = FALSE) {
    standardGeneric("get_net_cls_tab")
})
################################################################################

#' Get the cluster table 'cls_tab' from `norm` and `cls`.
#'
#' @param x An object of class `mina` with `norm` and `cls` defined.
#' @param uw By summing up the number of present components of each cluster
#' instead of relative abundances, default is FALSE.
#' @return The network cluster relative abundance table.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref", depth = 1000)
#' maize <- fit_tabs(maize)
#' maize <- adj(maize, method = "spearman")
#' maize <- net_cls(maize, method = "ap", cutoff = 0.5)
#' maize <- net_cls_tab(maize)
#' @export

setGeneric("net_cls_tab", function(x, uw = FALSE) {
    standardGeneric("net_cls_tab")
})

################################################################################

#' Inferring the network of different group of samples and test significance by
#' permutation.
#'
#' @param x An object of class `mina` with `norm` and `des` defined.
#' @param group The column name of descriptive file `des` for comparison.
#' @param ... Additional parameters.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- get_rep(maize, top = 5)
#' maize <- bs_pm(maize, group = "Compartment", per = 0.5)
#' @return The network bootstrap and permutation result.
#' @export

setGeneric("bs_pm", function(x, group, ...) {
    standardGeneric("bs_pm")
})

################################################################################

#' Calculate the network distance of `multi` and test the significance when
#' `perm` is defined.
#'
#' @param x An object of class `mina` with `multi` (and `perm` if `sig` is TRUE)
#' defined.
#' @param method The distance to be calculated, "spectra" and "Jaccard" are
#' available.
#' @param ... Additional parameters.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- get_rep(maize, top = 5)
#' maize <- bs_pm(maize, group = "Compartment")
#' maize <- net_dis(maize, method = "spectra", evk = 30)
#' @return The netowrk comparison result.
#' @export

setGeneric("net_dis", function(x, method, ...) {
    standardGeneric("net_dis")
})

################################################################################

#' Calculate the network distance of bootstrap and permutation when appliable.
#'
#' @param x The folder store the network inference results.
#' @param method The distance to be calculated, "spectra" and "Jaccard" are
#' available.
#' @param ... Additional parameters.
#' @return y The `mina` object with `dis_bs`, `dis_pm` and `dis_stat`.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- get_rep(maize, top = 5)
#' maize <- bs_pm(maize, group = "Compartment", individual = TRUE, out_dir =
#' "./individual_bs_pm/")
#' maize_stat1 <- net_dis_indi(x = "./individual_bs_pm/", method = "spectra")
#' maize_stat2 <- net_dis_indi(x = "./individual_bs_pm/", method = "Jaccard")
#' maize_stat3 <- net_dis_indi(x = "./individual_bs_pm/", method = "spectra",
#' evk = 100, skip = TRUE)
#' }
#' @export

setGeneric("net_dis_indi", function(x, method, ...) {
    standardGeneric("net_dis_indi")
})

##############################################################################

#' Get the representative community members by extracting the most abundant and
#' prevalent compositions.
#' @param x A quantitative matrix with samples in columns and compositions in
#' rows.
#'
#' @param top The percent of the most abundant and prevalent members.
#' @return The matrix with samples in columns and representative compositions in
#' rows.
#' @examples
#' data(maize_asv)
#' maize_asv_rep <- get_rep(maize_asv, top = 5)
#' @rdname get_rep-matrix
#' @export

setGeneric("get_rep", function(x, top = 5) {
    standardGeneric("get_rep")
})

###############################################################################

#' Visulization of network distance, average distances are used for tile plot.
#'
#' @import ggplot2
#' @param x An object of `mina` with slot `dis_stat` defined.
#' @param sig If `TRUE`, indicating significant distance with gold guild.
#' @param d The distance to be plotted, could be "BS" or "PM".
#' @return p The plotted figure.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment")
#' maize <- net_dis(maize, method = "Jaccard")
#' p <- net_dis_plot(maize)
#' @rdname net_dis_plot
#' @export

setGeneric("net_dis_plot", function(x, d = "BS", ...) {
    standardGeneric("net_dis_plot")
})

###############################################################################

#' Visulization of spectra network distance as PCoA.
#'
#' @import ggplot2
#' @import stringr
#' @param x The folder with all egv files generated by net_dis_indi().
#' @return p The plotted figure.
#' @examples
#' \dontrun{
#' data(maize)
#' norm(maize) <- maize_asv2
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment", individual = TRUE, out_dir =
#' "./individual_bs_pm/")
#' maize <- net_dis_indi("./individual_bs_pm/", method = "spectra", egv = TRUE,
#' dir = "./egv_folder/")
#' p <- net_dis_pcoa("./egv_folder/")
#' }
#' @rdname net_dis_pcoa
#' @export

setGeneric("net_dis_pcoa", function(x) {
    standardGeneric("net_dis_pcoa")
})

###############################################################################

#' Compare the node features between networks.
#'
#' @param x The folder with all network inference results generated by bs_pm()
#' @param cmp The compared feature of node, default `contrast`.
#' @param dir The directory to store the alculated node features.
#' @examples
#' \dontrun{
#' net_node_cmp("./individual_bs_pm/", f = "contrast", dir = "./")
#'}
#' @rdname net_node_cmp
#' @export

setGeneric("net_node_cmp", function(x, cmp = "contrast", dir = "./") {
    standardGeneric("net_node_cmp")
})

###############################################################################

#' Compare the group features between networks.
#'
#' @param x The folder with all network inference results generated by bs_pm()
#' @param cmp The compared feature of grp, default `contrast`.
#' @param dir The directory to store the alculated node features.
#' @param grp The table with group information.
#' @examples
#' \dontrun{
#' net_node_cmp("./individual_bs_pm/", f = "contrast", dir = "./", grp =
#' cls_tab(maize))
#'}
#' @rdname net_grp_cmp
#' @export

setGeneric("net_grp_cmp", function(x, cmp = "contrast", dir = "./", grp) {
    standardGeneric("net_grp_cmp")
})
