###############################################################################

#' Calculate the adjacacency matrix of @norm by correlation with `mina` class
#' object as input.
#'
#' @include all_classes.R all_generics.R
#' @import foreach, bigmemory, doMC
#' @param x An object of the class mina with @norm defined.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param threads (optional) The number of threads used for parallel running, 80
#' by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' x <- adj(x, method = "pearson")
#' x <- adj(x, method = "spearman", thread = 80, nblocks = 400)
#' x <- adj(x, method = "sparcc", thread = 40, nblocks = 200)
#' @return x An object of the class mina with @adj added.
#' @exportMethod adj

setMethod("adj", signature("mina", "ANY", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              stop("Must specify a `method`. See `? adj_method_list`")
          }
)

setMethod("adj", signature("mina", "character", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              x@adj <- adj(x@norm, method = method, threads = threads,
                           nblocks = nblocks)
              return(x)
          }
)

#' Calculate the adjacacency matrix of @norm by correlation with matrix as
#' input.
#'
#' @include all_classes.R all_generics.R
#' @import foreach, bigmemory, doMC
#' @param x An matrix for correlation calculation.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param threads (optional) The number of threads used for parallel running, 80
#' by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' y <- adj(x, method = "pearson", threads = 80, nblocks = 400)
#' @return y The adjacacency matrix.
#' @exportMethod adj

setMethod("adj", signature("matrix", "character", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              x <- t(x)

              if (method == "pearson" || method == "spearman") {
                  adj <- cor_par(x, method = method, threads = threads,
                                 nblocks = nblocks)
              }

              if (method == "sparcc") {
                  adj <- sparcc(x, threads = threads, nblocks = nblocks)
              }

              return(adj)
          }
)

#' Function for `pearson` and `spearman` correlation calculation. Modified from
#' https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
#' master/functions.community_similarity.R
#'
#' @include all_classes.R all_generics.R
#' @import foreach bigmemory doMC plyr
#' @importFrom parallel mclapply
#' @param x An matrix for correlation calculation.
#' @param method The correlation coeffient used for adjacacency matrix.
#' @param threads (optional) The number of threads used for parallel running, 80
#' by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' y <- cor_par(x, method = "pearson", threads = 80, nblocks = 400)
#' @return y The adjacacency matrix.
#' @keywords internal

cor_par <- function(x, method = "pearson", threads = 80, nblocks = 400) {
    use <- "na.or.complete"

    # Register cluster
    registerDoMC(cores = threads)

    # Preallocate blocks for parallel processing
    # => based on https://gist.github.com/bobthecat/5024079
    n_el <- ncol(x)

    if (n_el < 2 * nblocks) {
        x_rho <- cor(x, method = method, use = use)
        return(as.matrix(x_rho))
    }

    size_split <- floor(n_el / nblocks)

    if (size_split < 1) { size.split <- 1 }
    my_split <- list()
    length(my_split) <- nblocks
    my_split[1 : (nblocks - 1)] <- split(1 : (size_split * (nblocks - 1)),
                                         rep(1 : (nblocks - 1), each = size_split))
    my_split[[nblocks]] <- (size_split * (nblocks - 1)) : n_el

    dat_split <- mclapply(my_split, function(g) {
                                    as.matrix(x[, g])
                                    }, mc.cores = threads)

    # Get combinations of splits
    my_combs <- expand.grid(1 : length(my_split), 1:length(my_split))
    my_combs <- t(apply(my_combs, 1, sort))
    my_combs <- unique(my_combs)

    # Preallocate correlation matrix as big.matrix ("shared" in memory, so
    # accessible from w/in foreach loop)
    x_rho <- big.matrix(nrow = n_el, ncol = n_el, dimnames = list(colnames(x),
                                                                  colnames(x)),
                                                                  shared=T)
    x_rho_desc <- describe(x_rho)

    # Compute correlation matrix iterate through each block combination,
    # calculate matrix between blocks and store them in the preallocated matrix
    # on both symmetric sides of the diagonal.
    results <- foreach(i = 1 : nrow(my_combs)) %dopar% {
        # Get current combination and data
        curr_comb <- my_combs[i, ]
        g_1 <- my_split[[curr_comb[1]]]
        g_2 <- my_split[[curr_comb[2]]]
        curr_rho <- cor(x = dat_split[[curr_comb[1]]],
                        y = dat_split[[curr_comb[2]]], method=method)
        # Store
        curr_x_rho <- attach.big.matrix(x_rho_desc)
        curr_x_rho[g_1, g_2] <- curr_rho
        curr_x_rho[g_2, g_1] <- t(curr_rho)
        # Return
        TRUE
    }
    file.remove(list.files("/dev/shm/", full.name = T))
    return(as.matrix(x_rho))
}
