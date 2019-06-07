###############################################################################

#' Calculate the community dissimilarity / distance matrix of the input matrix.
#'
#' @include all_classes.R all_generics.R
#' @param x A matrix of the quantitative table.
#' @param method The dissimilarity / distance method used.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' y <- com_dis(x, method = "bray")
#' y <- com_dis(x, method = "tina", threads = 80, nblocks = 400)
#' @return y The dissimilarity / distance matrix.
#' @exportMethod com_dis

setMethod("com_dis", signature("matrix", "ANY", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              stop("Must specify a `method`, see `? com_dis_list`.")
          }
)

setMethod("com_dis", signature("matrix", "character", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              if (method == "tina"){
                  y <- tina(x, cor_method = "spearman", sim_method = "w_ja",
                            threads = threads, nblocks = nblocks)
              } else {
                  x <- t(x)
                  y <- dis_par(x, method = method, threads = threads,
                                   nblocks = nblocks)
              }
              return(y)
          }
)

###############################################################################

#' Calculate the community dissimilarity / distance matrix of @norm with `mina`
#' class object as input.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of the class mina with @norm defined.
#' @param method The dissimilarity / distance method used.
#' @param threads (optional) The number of threads used for parallel running.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix.
#' @examples
#' x <- com_dis(x, method = "bray")
#' x <- com_dis(x, method = "tina", threads = 40, nblocks = 200)
#' @return x The same `mina` object with @dis added.
#' @exportMethod com_dis

setMethod("com_dis", signature("mina", "ANY", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              stop("Must specify a `method`, see `? com_dis_list`.")
          }
)

setMethod("com_dis", signature("mina", "character", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              x@dis <- com_dis(x@norm, method = method, threads = threads,
                               nblocks = nblocks)
              return(x)
          }
)

###############################################################################

#' Function for `dis_par` dissimilarity / distance calculation. Modified from
#' https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
#' master/functions.community_similarity.R
#'
#' @include all_classes.R all_generics.R
#' @import foreach bigmemory doMC dplyr vegan
#' @importFrom parallel mclapply
#' @param x An matrix for dissimilarity / distance calculation.
#' @param method The method for dissimilarity / distance calculation.
#' @param threads (optional) The number of threads used for parallel running,
#' 80 by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' y <- dis_par(x, method = "bray", threads = 80, nblocks = 400)
#' @return y The dissimilarity / distance matrix.
#' @export

dis_par <- function(x, method = "bray", threads = 80, nblocks = 400) {
    use <- "na.or.complete"

    # Register cluster
    registerDoMC(cores = threads)

    nr <- nrow(x)
    if (nr < 2 * nblocks){
        return(as.matrix(vegdist(x, method = method, use = use)))
    }

    size_split <- floor(nr / nblocks)
    size_split <- size_split + 1
    nblocks <- floor(nr / size_split)
    nblocks <- nblocks + 1

    my_split <- list()
    length(my_split) <- nblocks
    my_split[1 : (nblocks - 1)] <- split(1 : (size_split * (nblocks - 1)),
                                        rep(1 : (nblocks - 1), each = size_split))
    my_split[[nblocks]] <- (size_split * (nblocks - 1)) : nr

    dat_split <- mclapply(my_split, function(g) { as.matrix(x[, g]) },
                          mc.cores = threads)

    # Get combinations of splits
    my_combs <- expand.grid(1 : length(my_split), 1:length(my_split))
    my_combs <- t(apply(my_combs, 1, sort))
    my_combs <- unique(my_combs)

    # Preallocate dissimilarity / distance matrix as big.matrix
    x_dis <- big.matrix(nrow = nr, ncol = nr,
                        dimnames = list(rownames(x), rownames(x)), shared=T)
    x_dis_desc <- describe(x_dis)

    # Compute dissimilarity / distance matrix iterate through each block
    # combination, calculate matrix between blocks and store them in the
    # preallocated matrix on bxh symmetric sides of the diagonal.

    results <- foreach(i = 1 : nrow(my_combs)) %dopar% {
        # Get current combination and data
        curr_comb <- my_combs[i, ]

        g_1 <- my_split[[curr_comb[1]]]
        g_2 <- my_split[[curr_comb[2]]]
        data_1 <- dat_split[[curr_comb[1]]]
        data_2 <- dat_split[[curr_comb[2]]]

        if (curr_comb[1] == curr_comb[2]) {
            x_dis <- as.matrix(vegdist(data_1, method = method))
        } else {
            x_dis <- as.matrix(vegdist(rbind(data_1, data_2), method = method))
            n1 <- nrow(data_1)
            n2 <- nrow(dis)
            x_dis <- x_dis[1 : n1, (n1 + 1) : n2]
        }
        # Store
        curr_x_dis <- attach.big.matrix(x_dis_desc)
        curr_x_dis[g_1, g_2] <- curr_dis
        curr_x_dis[g_2, g_1] <- t(curr_dis)
        # Return
        TRUE
    }

    file.remove(list.files("/dev/shm/", full.name = T))
    return(as.matrix(x_dis))
}


###############################################################################

#' Function for `tina` dissimilarity / distance calculation. Modified from
#' https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
#' master/functions.community_similarity.R
#' Pearson / Spearman could be used for calculating correlation and weighted /
#' unweighted Jaccard could be used for the calculation of similarity.
#'
#' @include all_classes.R all_generics.R
#' @import foreach bigmemory doMC dplyr vegan
#' @importFrom parallel mclapply
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

tina <- function(x, cor_method = "spearman", sim_method = "w_ja",
                 threads = 80, nblocks = 400) {
    x_sparcc <- sparcc(x, threads = threads, nblocks = nblocks)
    tmp.S <- cor_par(x_sparcc, method = cor_method, threads = threads,
                     nblocks = nblocks)
    Cij <- 0.5 * (tmp.S + 1)
    t <- sim_par(x, Cij, sim_method = sim_method, threads = threads,
                 nblocks = nblocks)
    t[t < 0] <- 0
    return(t)
}

###############################################################################

#' Function for community similarity calculation used by `tina`, modified from
#' https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
#' master/functions.community_similarity.R
#'
#' @include all_classes.R all_generics.R
#' @import plyr foreach bigmemory doMC
#' @importFrom parallel mclapply
#' @param x An quantitative matrix.
#' @param y The Cij matrix, which is correlation matrix of adjusted sparcc
#' matrix of x.
#' @param sim_method The method for similarity, "w_ja" and "uw_ja" are
#' available for weighted and unweighted Jaccard similarity respectively.
#' @param threads (optional) The number of threads used for parallel running,
#' 80 by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' s <- sim_par(x, y, sim_method = "w_ja", threads = 80, nblocks = 400)
#' @return s The output similarity matrix.
#' @keywords internal

sim_par <- function(x, y, sim_method = "w_ja", threads = 80, nblocks = 400) {
    registerDoMC(cores = threads)

    samples <- colnames(x)
    x_occ_list <- apply((x > 0), 2, which)
    x_count_list <- alply(x, .margins=2,
                           .fun=function(x_vec) { x_vec[x_vec > 0] },
                           .parallel=T)
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

    file.remove(list.files("/dev/shm/", full.name=T))

    s <- do.call("rbind", cs_list)
    rownames(s) <- colnames(s) <- samples
    return(s)
}


###############################################################################

#' List of dissimilarity / distance supported in \code{\link[mina]{com_dis}}
#'
#' Dissimilarity / distance should be specified by exact string match.
#'
#' @format A list of character vectors indicate the dissimilarity / distance
#' method used.
#'
#' \describe{
#'   \item{\code{tina}}{TINA from Schmidt_et_al_2016}
#'
#'   \item{\code{manhattan}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{euclidean}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{canberra}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{bray}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{kulczynski}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{jaccard}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{gower}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{altGower}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{morisita}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{horn}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{mountford}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{raup}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{binomial}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{chao}}{ from \code{\link[vegan]{vegdist}}}
#'   \item{\code{cao}}{ from \code{\link[vegan]{vegdist}}}
#' }
#'
#' @export
#' @examples
#' com_dis_list

com_dis_list <- list(
    # The methods supported by vegan::vegdist function.
    vegdist    = c("manhattan", "euclidean", "canberra", "bray",
                   "kulczynski", "jaccard", "gower", "altGower", "morisita",
                   "horn", "mountford", "raup" , "binomial", "chao", "cao"),
    TINA       = "tina",
    design_dis = "ANY"
)
