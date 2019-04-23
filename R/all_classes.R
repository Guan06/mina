#' Class "mina" includes the quantitative table and descriptive table.
#'
#' @name mina-class
#' @aliases mina-class show, mina-method comDis
#' @docType class
#' @slot tab the quantitative table of the dataset.
#' @slot des_tab the descriptive table of the samples listed in Tab.
#' @slot norm the normalized quantitative table of Tab.
#' @slot dis the distance / dissimilarity matrix between samples in Tab and desTab.
#' @slot dmr the list of dimensionality reduction result, includes points and variance.
#' @slot adj the adjacency matrix between pairwise compositions (e.g. OTUs/ASVs)
#' @slot cls the cluster information for each composition
#' @slot cls_tab the cluster quantitative table
#' @slot dis_net the distance between networks of different environmental communities.
#' @slot multi the list of subsampled adjacency matrices for each environment.
#' @slot perm the list of permutated adjacency matrices for each pairwise environmental comparison.
#' @slot sig the average distance between subsampled environmental community networks and corresponding
#'          significance.
#'
#' @examples
#' new(mina, tab = "otu_table.txt", des_tab = "design.txt")
#' @author Rui Guan \url{https://github.com/Guan06}
#' @export
setClass('mina',
         representation(tab = "matrix",
                   des_tab = "data.frame",
                   norm = "matrix",
                   dis = "matrix",
                   dmr = "list",
                   adj = "matrix",
                   cls = "data.frame",
                   cls_tab = "matrix",
                   dis_net = "data.frame",
                   multi = "list",
                   perm = "list",
                   sig = "data.frame"
         ),
         prototype(tab = NULL, des_tab = NULL, norm = NULL,
                   dis = NULL, dmr = NULL, adj = NULL,
                   cls = NULL, cls_tab = NULL, dis_net = NULL,
                   multi = NULL, perm = NULL, sig = NULL)
)
