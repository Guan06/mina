#' Class "mina" includes the quantitative table and descriptive table.
#'
#' @name mina-class
#' @aliases mina-class show, mina-method comDis
#' @docType class
#' @slot Tab the quantitative table of the dataset.
#' @slot desTab the descriptive table of the samples listed in Tab.
#' @slot normTab the normalized quantitative table of Tab.
#' @slot dis the distance / dissimilarity matrix between samples in Tab and desTab.
#' @slot dmr the list of dimensionality reduction result, includes points and variance.
#' @slot adj the adjacency matrix between pairwise compositions (e.g. OTUs/ASVs)
#' @slot cls the cluster information for each composition
#' @slot clsTab the cluster quantitative table
#' @slot disNet the distance between networks of different environmental communities.
#' @slot multi the list of subsampled adjacency matrices for each environment.
#' @slot perm the list of permutated adjacency matrices for each pairwise environmental comparison.
#' @slot sig the average distance between subsampled environmental community networks and corresponding
#'          significance.
#'
#' @examples
#' new(mina, Tab = "otu_table.txt", desTab = "design.txt")
#' @author Rui Guan \url{https://github.com/Guan06}
#' @export
setClass('mina',
         representation(Tab = "matrix",
                   desTab = "data.frame",
                   normTab = "matrix",
                   dis = "matrix",
                   dmr = "list",
                   adj = "matrix",
                   cls = "data.frame",
                   clsTab = "matrix",
                   disNet = "data.frame",
                   multi = "list",
                   perm = "list",
                   sig = "data.frame"
         ),
         prototype(Tab = NULL, desTab = NULL, normTab = NULL,
                   dis = NULL, dmr = NULL, adj = NULL,
                   cls = NULL, clsTab = NULL, disNet = NULL,
                   multi = NULL, perm = NULL, sig = NULL)
)
