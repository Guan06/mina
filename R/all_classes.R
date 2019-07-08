###############################################################################

#' Modified from https://github.com/joey711/phyloseq/blob/master/R/allClasses.R
#' Use setClassUnion to define the unholy NULL-data union as a virtual class.
#' This is a way of dealing with the expected scenarios in which one or more of
#' the component data classes is not available, in which case NULL will be used
#' instead.

#' @keywords internal

setClassUnion("mat_or_NULL", c("matrix", "NULL"))

#' @keywords internal

setClassUnion("df_or_NULL", c("data.frame", "NULL"))

#' @keywords internal

setClassUnion("lst_or_NULL", c("list", "NULL"))

###############################################################################

#' Class "mina" includes the quantitative table and descriptive table.
#'
#' @name mina-class
#' @aliases mina-class
#' @docType class
#' @slot tab The quantitative table of the dataset.
#' @slot des The descriptive table of the samples listed in @tab.
#' @slot norm The normalized quantitative table of @tab.
#' @slot dis The distance / dissimilarity matrix between samples in @tab.
#' @slot dmr The list of dimensionality reduction result, includes points and
#' variance.
#' @slot adj The adjacency matrix between pairwise compositions (e.g. OTUs/ASVs)
#' @slot cls The cluster information for each composition.
#' @slot cls_tab The cluster quantitative table.
#'
#' @slot multi The list of subsampled adjacency matrices for each environment.
#' @slot perm The list of permutated adjacency matrices for each pairwise
#' environmental comparison.
#' @slot dis_bs The distance between networks of different environmental
#' communities.
#' @slot dis_pm The distance between networks of permutated groups.
#' @slot dis_stat The average distance between subsampled environmental community
#' networks, permutated networks and corresponding significance.
#'
#' @examples
#' new("mina", tab = as.matrix(maize_asv), des = maize_des)
#' @author Rui Guan \url{https://github.com/Guan06}
#' @exportClass mina

setClass("mina",
         representation(tab = "mat_or_NULL",
                   des = "df_or_NULL",
                   norm = "mat_or_NULL",
                   dis = "mat_or_NULL",
                   dmr = "lst_or_NULL",
                   adj = "mat_or_NULL",
                   cls = "df_or_NULL",
                   cls_tab = "mat_or_NULL",

                   multi = "lst_or_NULL",
                   perm = "lst_or_NULL",
                   dis_bs = "df_or_NULL",
                   dis_pm = "df_or_NULL",
                   dis_stat = "df_or_NULL"
         ),
         prototype(tab = NULL, des = NULL, norm = NULL,
                   dis = NULL, dmr = NULL, adj = NULL,
                   cls = NULL, cls_tab = NULL,
                   multi = NULL, perm = NULL, dis_bs = NULL, dis_pm = NULL,
                   dis_stat = NULL)
)

################################################################################
