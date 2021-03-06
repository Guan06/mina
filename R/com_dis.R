###############################################################################

#' Calculate the community dissimilarity / distance matrix of the input matrix.
#'
#' @include all_classes.R all_generics.R
#' @importFrom parallelDist parDist
#' @param x A matrix of the quantitative table.
#' @param method The dissimilarity / distance method used, default `bray`.
#' @param threads (optional, only needed when method == "tina") The number of
#' threads used for parallel running.
#' @param nblocks (optional, only needed when method == "tina") The number of
#' row / column for splitted sub-matrix.
#' @param ... Additional parameters.
#' @examples
#' asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
#' asv_dis <- com_dis(asv_norm, method = "bray")
#' @return y The dissimilarity / distance matrix.
#' @rdname com_dis-matrix
#' @exportMethod com_dis

setMethod("com_dis", signature("matrix", "ANY"),
          function(x, method = "bray", threads = 80, nblocks = 400, ...) {
              stop("Must specify a `method`, see `? com_dis_list`.")
          }
)

###############################################################################

#' @rdname com_dis-matrix
#' @exportMethod com_dis

setMethod("com_dis", signature("matrix", "character"),
          function(x, method = "bray", threads = 80, nblocks = 400, ...) {
              stopifnot(
                        method %in% unlist(com_dis_list),
                        is.numeric(c(threads, nblocks))
              )
              if (method == "tina"){
                  y <- tina(x, cor_method = "spearman", sim_method = "w_ja",
                            threads = threads, nblocks = nblocks)
              } else if (method == "Jaccard") {
                  x <- t(x)
                  y <- as.matrix(parDist(x, method = "bray"))
                  y <- 2 * y / (1 + y)
              }else {
                  x <- t(x)
                  y <- as.matrix(parDist(x, method = method))
              }
              return(y)
          }
)

###############################################################################

#' Calculate the community dissimilarity / distance matrix of `norm` with `mina`
#' class object as input.
#'
#' @include all_classes.R all_generics.R
#' @importFrom parallelDist parDist
#' @param x An object of the class `mina` with `norm` defined.
#' @param method The dissimilarity / distance method used, default `bray`.
#' @param threads (optional, only needed when method == "tina") The number of
#' threads used for parallel running.
#' @param nblocks (optional, only needed when method == "tina") The number of
#' row / column for splitted sub-matrix.
#' @param ... Additional parameters.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "total")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' @return x The same `mina` object with @dis added.
#' @rdname com_dis-mina
#' @exportMethod com_dis

setMethod("com_dis", signature("mina", "ANY"),
          function(x, method = "bray", threads = 80, nblocks = 400, ...) {
              stop("Must specify a `method`, see `? com_dis_list`.")
          }
)

###############################################################################

#' @rdname com_dis-mina
#' @exportMethod com_dis
setMethod("com_dis", signature("mina", "character"),
          function(x, method = "bray", threads = 80, nblocks = 400, ...) {
              stopifnot(
                        method %in% unlist(com_dis_list),
                        is.numeric(c(threads, nblocks))
              )
              dis(x) <- com_dis(norm(x), method = method, threads = threads,
                               nblocks = nblocks)
              return(x)
          }
)

###############################################################################

#' Function for `tina` dissimilarity calculation. Modified from Schmidt et al.,
#' 2016. Person and Spearman could be used for correlation and weighted and
#' unweighted Jaccard could be used for similarity calculation.
#'
#' @include all_classes.R all_generics.R
#' @param x A matrix for dissimilarity calculation.
#' @param cor_method The method for correlation, "pearson" and "spearman" are
#' available.
#' @param sim_method The method for similarity, "w_ja" and "uw_ja" are
#' available for weighted and unweighted Jaccard similarity respectively.
#' @param threads The number of threads used for parallel running, 80 by
#' default.
#' @param nblocks The number of row and column for splitted sub-matrix, 400 by
#' default.
#' @param ... Additional parameters.
#' @examples
#' \dontrun{
#' asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000)
#' asv_dis <- com_dis(asv_norm, method = "bray")
#' asv_dis <- com_dis(asv_norm, method = "tina", threads = 8, nblocks = 40)
#' asv_tina <- tina(asv_norm, cor_method = "spearman", sim_method = "w_ja",
#' threads = 8, nblocks = 40)
#' }
#' @return t The output `tina` dissimilarity matrix.
#' @exportMethod tina

setMethod("tina", signature("matrix"),
    function(x, cor_method = "spearman", sim_method = "w_ja",
                 threads = 80, nblocks = 400, ...) {
        stopifnot(
                cor_method %in% c("spearman", "pearson", "sparcc"),
                sim_method %in% c("w_ja", "uw_ja"),
                is.numeric(c(threads, nblocks))
        )
    x_sparcc <- sparcc(x, threads = threads, nblocks = nblocks)
    tmp.S <- adj(x_sparcc, method = cor_method)

    Cij <- 0.5 * (tmp.S + 1)
    t <- sim_par(x, Cij, sim_method = sim_method, threads = threads,
                 nblocks = nblocks)
    t[t < 0] <- 0
    return(t)
})

###############################################################################

#' Function for community similarity calculation used by `tina`, modified from
#' https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
#' master/functions.community_similarity.R
#'
#' @include all_classes.R all_generics.R
#' @importFrom parallel mclapply
#' @importFrom foreach foreach
#' @importFrom bigmemory big.matrix
#' @importFrom bigmemory attach.big.matrix
#' @importFrom plyr alply
#'
#' @param x An quantitative matrix.
#' @param y The Cij matrix, which is correlation matrix of adjusted sparcc
#' matrix of x.
#' @param sim_method The method for similarity, "w_ja" and "uw_ja" are
#' available for weighted and unweighted Jaccard similarity respectively.
#' @param threads The number of threads used for parallel running, 80 by
#' default.
#' @param nblocks The number of row / column for splitted sub-matrix, 400 by
#' default.
#' @examples
#' \dontrun{
#' data(maize_asv)
#' maize_tab <- maize_asv[1 : 1000, 1 : 200]
#' asv <- norm_tab(maize_tab, method = "raref", depth = 100)
#' asv[is.na(asv)] <- 0
#' asv_sparcc <- sparcc(asv, threads = 8, nblocks = 40)
#' tmp.S <- adj(asv_sparcc, method = "spearman")
#' y <- 0.5 * (tmp.S + 1)
#' s <- sim_par(asv_sparcc, y, sim_method = "w_ja", threads = 8, nblocks = 40)
#' }
#' @return s The output similarity matrix.
#' @keywords internal

sim_par <- function(x, y, sim_method = "w_ja", threads = 80, nblocks = 400) {
    doMC::registerDoMC(cores = threads)

    samples <- colnames(x)
    x_occ_list <- apply((x > 0), 2, which)
    x_count_list <- alply(x, .margins=2,
                           .fun=function(x_vec) { x_vec[x_vec > 0] },
                           .parallel=TRUE)
    names(x_occ_list) <- names(x_count_list) <- samples

    if (sim_method == "uw_ja") {
        smpl_csums <- mclapply(x_occ_list,
                                 function(a_l) {
                                     colSums(y[a_l, ]) / length(a_l)
                                 }, mc.cores = threads)
        names(smpl_csums) <- samples

        smpl_self <- lapply(samples,
                            function(a) {
                                al <- x_occ_list[[a]]
                                sum(smpl_csums[[a]][al]) / length(al)
                            } )
        smpl_self <- unlist(smpl_self)
        names(smpl_self) <- samples

        cs_list <- mclapply(samples,
                            function(a) {
                                a_sums <- smpl_csums[[a]]
                                a_self <- smpl_self[a]
                                unlist(lapply(samples,
                                              function(b) {
                                                  b_list <- x_occ_list[[b]]
                                                  b_self <- smpl_self[b]
                                                  m <- sum(a_sums[b_list])
                                                  ab <- sqrt(a_self * b_self)
                                                  n <- length(b_list) * ab
                                                  1 - m / n
                                              }
                                             )
                                )
                            }, mc.cores = threads)
    }


    if (sim_method == "w_ja") {
        x_rel_count <- lapply(x_count_list,
                              function(a_count) { a_count / sum(a_count) } )

        smpl_self <- unlist(mclapply(x_rel_count,
                                     function(a_rel) {
                                         a_list <- names(a_rel)
                                         a_y <- y[a_list, a_list]
                                         sum(a_y * outer(a_rel, a_rel))
                                     }, mc.cores = threads))

        cs_list <- mclapply(samples,
                            function(a) {
                                a_rel <- x_rel_count[[a]]
                                a_list <- names(a_rel)
                                a_y <- a_rel * y[a_list,]
                                unlist(lapply(samples,
                                              function(b) {
                                                  b_rel <- x_rel_count[[b]]
                                                  b_list <- names(b_rel)
                                                  curr_y <- a_y[, b_list]
                                                  m <- sum(b_rel * t(curr_y))
                                                  a <- smpl_self[a]
                                                  b <- smpl_self[b]
                                                  n <- sqrt(a * b)
                                                  1 - m / n
                                              } ) )
                            }, mc.cores = threads)
    }

    file.remove(list.files("/dev/shm/", full.names = TRUE))

    s <- do.call("rbind", cs_list)
    rownames(s) <- colnames(s) <- samples
    return(s)
}

###############################################################################


###############################################################################

#' List of dissimilarity / distance supported in \code{\link[mina]{com_dis}}.
#' Dissimilarity / distance should be specified by exact string match.
#'
#' @format A list of character vectors indicate the dissimilarity / distance
#' method used.
#'
#' \describe{
#'   \item{\code{tina}}{ TINA from Schmidt_et_al_2016 }
#'
#'   \item{\code{Jaccard}}{ Jaccard defined by \code{\link[vegan]{vegdist}} }
#'
#'   \item{weighted}{ Dissimilarity / distance method for weighted matrix: }
#'   \item{\code{bhjattacharyya}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{canberra}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{bray}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{chord}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{divergence}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{euclidean}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{fJaccard}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{geodesic}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{hellinger}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{kullback}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{manhattan}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{maximum}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{minkowski}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{podani}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{soergel}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{wave}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{whittaker}}{ from \code{\link[parallelDist]{parDist}} }
#'
#'   \item{unweighted}{ Dissimilarity / Distance for unweighted matrix: }
#'   \item{\code{binary}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{braun-blanquet}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{consine}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{dice}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{fager}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{faith}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{hamman}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{hamming}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{kulczynski1}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{kulczynski2}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{michael}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{mountford}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{mozley}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{ochiai}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{phi}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{russel}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{simple matching}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{simpson}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{stiles}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{tanimoto}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{yule}}{ from \code{\link[parallelDist]{parDist}} }
#'   \item{\code{yule2}}{ from \code{\link[parallelDist]{parDist}} }
#'
#' }
#'
#' @export
#' @examples
#' ? com_dis_list

com_dis_list <- list(
    # Dissimilarity / distance implemented in parallelDist
    par_dist    = c("bhjattacharyya", "bray", "canberra", "chord", "divergence",
                   "euclidean", "fJaccard", "geodesic", "hellinger", "kullback",
                   "manhattan", "maximum", "minkowski", "podani", "soergel",
                   "wave", "whittaker"),

    # Dissimilarity / distance implemented in parallelDist for binary matrix
    par_dist_bi = c("binary", "braun-blanquet", "consine", "dice", "fager",
                    "faith", "hamman", "hamming", "kulczynski1", "kulczynski2",
                    "michael", "mountford", "mozley", "ochiai", "phi", "russel",
                    "simple matching", "simpson", "stiles", "tanimoto", "yule",
                    "yule2"),
    TINA        = "tina",
    Jaccard     = "Jaccard"
)
