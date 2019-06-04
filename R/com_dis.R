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
#' y <- com_dis(x, method = "tina", threads = 40, nblocks = 200)
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
                  dis <- tina(x, threads = threads, nblocks = nblocks)
              } else {
                  x <- t(x)
                  dis <- dis_par(x, method = method, threads = threads,
                                   nblocks = nbolcks)
              }
              return(dis)
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

#' Function for dis_par dissimilarity / distance calculation. Modified from
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
#' @keywords internal

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
    # preallocated matrix on both symmetric sides of the diagonal.

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

#' List of dissimilarity / distance supported in \code{\link[mina]{com_dis}}
#'
#' Dissimilarity / distance should be specified by exact string match.
#'
#' @format A list of character vectors.
#'
#' \describe{
#'   \item{\code{tina}}{TINA from Schmidt_et_al_2016}
#'
#'   \item{\code{manhattan}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{euclidean}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{canberra}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{bray}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{kulczynski}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{jaccard}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{gower}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{altGower}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{morisita}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{horn}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{mountford}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{raup}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{binomial}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{chao}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{cao}}{\code{\link[vegan]{vegdist}}}
#' }
#'
#' @export
#'
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
