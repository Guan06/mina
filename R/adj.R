###############################################################################

#' Calculate the adjacacency matrix of @norm by correlation with `mina` class
#' object as input.
#'
#' @include all_classes.R all_generics.R
#' @import foreach bigmemory doMC
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
#' @return x The same `mina` object with @adj added.
#' @exportMethod adj

setMethod("adj", signature("mina", "ANY", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              stop("Must specify a `method`, see `? adj_method_list`.")
          }
)

setMethod("adj", signature("mina", "character", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              x@adj <- adj(x@norm, method = method, threads = threads,
                           nblocks = nblocks)
              return(x)
          }
)

###############################################################################

#' Calculate the adjacacency matrix of @norm by correlation with matrix as
#' input.
#'
#' @include all_classes.R all_generics.R
#' @import foreach bigmemory doMC
#' @param x An matrix for correlation / adjacency matrix calculation.
#' @param method The correlation coefficient used for adjacacency matrix.
#' @param threads (optional) The number of threads used for parallel running, 80
#' by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' y <- adj(x, method = "pearson", threads = 80, nblocks = 400)
#' @return y The adjacacency matrix.
#' @exportMethod adj

setMethod("adj", signature("matrix", "ANY", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              stop("Must specify a `method`, see `? adj_method_list`.")
          }
)

setMethod("adj", signature("matrix", "character", "ANY", "ANY"),
          function(x, method, threads = 80, nblocks = 400) {
              if (method == "pearson" || method == "spearman") {
                  x <- t(x)
                  adj <- cor_par(x, method = method, threads = threads,
                                 nblocks = nblocks)
              }

              if (method == "sparcc") {
                  adj <- sparcc(x, threads = threads, nblocks = nblocks)
              }

              return(adj)
          }
)

###############################################################################

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
        return(as.matrix(cor(x, method = method, use = use)))
    }

    size_split <- floor(n_el / nblocks)

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

    # Preallocate correlation matrix as big.matrix
    x_cor <- big.matrix(nrow = n_el, ncol = n_el,
                        dimnames = list(colnames(x), colnames(x)), shared=T)
    x_cor_desc <- describe(x_cor)

    # Compute correlation matrix iterate through each block combination,
    # calculate matrix between blocks and store them in the preallocated matrix
    # on both symmetric sides of the diagonal.
    results <- foreach(i = 1 : nrow(my_combs)) %dopar% {
        # Get current combination and data
        curr_comb <- my_combs[i, ]
        g_1 <- my_split[[curr_comb[1]]]
        g_2 <- my_split[[curr_comb[2]]]
        curr_cor <- cor(x = dat_split[[curr_comb[1]]],
                        y = dat_split[[curr_comb[2]]], method=method)
        # Store
        curr_x_cor <- attach.big.matrix(x_cor_desc)
        curr_x_cor[g_1, g_2] <- curr_cor
        curr_x_cor[g_2, g_1] <- t(curr_cor)
        # Return
        TRUE
    }
    file.remove(list.files("/dev/shm/", full.name = T))
    return(as.matrix(x_cor))
}

###############################################################################

#' Function for `sparcc` correlation calculation. Modified from
#' https://github.com/defleury/Schmidt_et_al_2016_community_similarity/blob/
#' master/functions.community_similarity.R
#'
#' @include all_classes.R all_generics.R
#' @import foreach bigmemory doMC plyr
#' @importFrom parallel mclapply
#' @param x An matrix for correlation calculation.
#' @param threads (optional) The number of threads used for parallel running, 80
#' by default.
#' @param nblocks (optional) The number of row / column for splitted sub-matrix,
#' 400 by default.
#' @examples
#' y <- sparcc(x, threads = 80, nblocks = 400)
#' @return y The adjacency matrix.
#' @keywords internal

sparcc <- function(x, threads = 80, nblocks = 400) {
    use <- "na.or.complete"

    # add pseudocount to avoid issues with 0 in log-space
    pseudocount <- 10^-6
    x <- x + pseudocount

    # Register cluster
    registerDoMC(cores = threads)

    # Preallocate blocks for parallel processing
    # => based on https://gist.github.com/bobthecat/5024079
    nr <- nrow(x)
    size_split <- floor(nr / nblocks)

    my_split <- list()
    length(my_split) <- nblocks
    my_split[1 : (nblocks - 1)] <- split(1 : (size_split * (nblocks - 1)),
                                rep(1 : (nblocks - 1), each = size_split))
    my_split[[nblocks]] <- (size_split * (nblocks - 1)) : nr

    dat_split <- mclapply(my_split, function(g) {
                                    as.matrix(x[, g])
                                    }, mc.cores = threads)

    # Get combinations of splits
    my_combs <- expand.grid(1 : length(my_split), 1:length(my_split))
    my_combs <- t(apply(my_combs, 1, sort))
    my_combs <- unique(my_combs)

    # Preallocate T matrix as big.matrix
    x_cor <- big.matrix(nrow = nr, ncol = nr,
                        dimnames = list(rownames(x), rownames(x)), shared=T)
    x_cor_desc <- describe(x_cor)

    # Compute Aitchinson's T matrix iterate through each block combination,
    # calculate matrix between blocks and store them in the preallocated matrix
    # on both symmetric sides of the diagonal.
    results <- foreach(i = 1 : nrow(my_combs)) %dopar% {
        # Get current combination and data
        curr_comb <- my_combs[i, ]
        g_1 <- my_split[[curr_comb[1]]]
        g_2 <- my_split[[curr_comb[2]]]
        data_1 <- dat_split[[curr_comb[1]]]
        data_2 <- dat_split[[curr_comb[2]]]

        curr_cor <- apply(data_1, 1, function(x) {
                        apply(data_2, 1, function(y) { var(log(x / y)) })
                    })
        # Store
        curr_x_cor <- attach.big.matrix(x_cor_desc)
        curr_x_cor[g_1, g_2] <- curr_cor
        curr_x_cor[g_2, g_1] <- t(curr_cor)
        # Return
        TRUE
    }

    # Compute component variations t_i
    var <- colsum(x_cor)
    # Estimate component variances ("omega_i") from t_i by solving a linear
    # equation system
    mat_a <- matrix(data = 1, nrow = nr, ncol = nr)
    diag(mat_a) <- nr - 1
    omega <- sqrt(solve(a = mat_a, b = var))
    spa <- foreach(i = 1 : nr, .combine = 'rbind', .multicombine = T) %dopar% {
              (omega[i] ^ 2 + omega ^ 2 - x_cor[i, ]) / (2 * omega[i] * omega)
    }

    rownames(spa) <- rownames(x)
    file.remove(list.files("/dev/shm/", full.name = T))
    return(spa)
}

###############################################################################

#' List of adjacacency matix calculation methods / correlations supported in
#' \code{\link[mina]{adj}}
#'
#' Correlation methods should be specified by exact string match.
#' @format A list of character vectors.
#' \describe{
#'    \item{pearson}{
#'        Pearson correlation.
#'    }
#'    \item{spearman}{
#'        Spearman correlation.
#'    }
#'    \item{sparcc}{
#'        SparCC correlation by spearman.
#'    }
#'
#' }
#' @seealso \code{\link[mina]{adj}}
#' @export
#' @examples
#' ? adj_method_list

adj_method_list <- list(
    pearson = "pearson",
    spearman = "spearman",
    sparcc = "sparcc"
)