###############################################################################
# Modified from https://github.com/joey711/phyloseq/blob/master/R/allClasses.R
# Use setClassUnion to define the unholy NULL-data union as a virtual class.
# This is a way of dealing with the expected scenarios in which one or more of
# the component data classes is not available, in which case NULL will be used
# instead.

#' @keywords internal
setClassUnion("mat_or_NULL", c("matrix", "NULL"))
#'
#' @keywords internal
setClassUnion("df_or_NULL", c("data.frame", "NULL"))
#'
#' @keywords internal
setClassUnion("lst_or_NULL", c("list", "NULL"))

###############################################################################

#' Class "mina" includes the quantitative table and descriptive table.
#'
#' @name mina-class
#' @aliases mina-class show, mina-method comDis
#' @docType class
#' @slot tab The quantitative table of the dataset.
#' @slot des_tab The descriptive table of the samples listed in @tab.
#' @slot norm The normalized quantitative table of @tab.
#' @slot dis The distance / dissimilarity matrix between samples in @tab.
#' @slot dmr The list of dimensionality reduction result, includes points and
#' variance.
#' @slot adj The adjacency matrix between pairwise compositions (e.g. OTUs/ASVs)
#' @slot cls The cluster information for each composition.
#' @slot cls_tab The cluster quantitative table.
#' @slot dis_net The distance between networks of different environmental
#' communities.
#' @slot multi The list of subsampled adjacency matrices for each environment.
#' @slot perm The list of permutated adjacency matrices for each pairwise
#' environmental comparison.
#' @slot sig The average distance between subsampled environmental community
#' networks and corresponding significance.
#'
#' @examples
#' new("mina", tab = as.matrix(maize_asv), des_tab = maize_des)
#' @author Rui Guan \url{https://github.com/Guan06}
#' @exportClass mina
setClass("mina",
         representation(tab = "mat_or_NULL",
                   des_tab = "df_or_NULL",
                   norm = "mat_or_NULL",
                   dis = "mat_or_NULL",
                   dmr = "lst_or_NULL",
                   adj = "mat_or_NULL",
                   cls = "df_or_NULL",
                   cls_tab = "mat_or_NULL",
                   dis_net = "df_or_NULL",
                   multi = "lst_or_NULL",
                   perm = "lst_or_NULL",
                   sig = "df_or_NULL"
         ),
         prototype(tab = NULL, des_tab = NULL, norm = NULL,
                   dis = NULL, dmr = NULL, adj = NULL,
                   cls = NULL, cls_tab = NULL, dis_net = NULL,
                   multi = NULL, perm = NULL, sig = NULL)
)

################################################################################
