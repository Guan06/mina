###############################################################################

#' Normalize the quantitative matrix.
#'
#' @include all_classes.R all_generics.R
#' @param x A quantitative matrix with samples in columns and compositions in
#' rows.
#' @param method The method used for normalization.
#' @param depth The depth for rarefying, 1000 by default.
#' @param replace Whether to sample with replacement (\code{TRUE} by default)
#' or without replacement (\code{FALSE}) when using method `raref`.
#' @param multi Rarefy the table for multiple times, FALSE by default, indicate
#' the times of rarefaction want to be repeated, only validate for rarefaction.
#' @examples
#' data(maize_asv2)
#' maize_asv_norm <- norm_tab(maize_asv2, method = "total")
#' maize_asv_norm <- norm_tab(maize_asv2, method = "raref", depth = 1000,
#' replace = TRUE, multi = 3)
#' @return x_norm Normalized matrix of the quantitative table.
#' @rdname norm_tab-matrix
#' @exportMethod norm_tab

setMethod("norm_tab", signature("matrix", "character", "ANY", "ANY", "ANY"),
          function(x, method, depth = 1000, replace = TRUE, multi = FALSE) {
              if (method == "raref") {
                  if (multi == FALSE) {
                    x_norm <- norm_by_raref(x, depth = depth, replace = replace)
                  } else {
                      multi_norm <- norm_by_raref(x, depth = depth,
                                                  replace = replace)
                      for ( i in 2 : multi ){
                          multi_norm <- multi_norm + norm_by_raref(x,
                                                        depth = depth,
                                                        replace = replace)
                      }
                      multi_norm <- multi_norm[rowSums(multi_norm) > 0, ]
                      x_norm <- multi_norm / multi
                  }
              }
              if (method == "total") {
                  x_norm <- norm_by_total(x)
              }
              x_norm <- x_norm[rowSums(x_norm) > 0, ]
              return(as.matrix(x_norm))
          }
)

################################################################################

#' Normalize the quantitative table with mina input.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of the class mina with @tab defined.
#' @param method The method used for normalization.
#' @param depth The depth for subsampling by rarefying, 1000 by default.
#' @param replace Whether to sample with replacement (\code{TRUE} by default) or
#' without replacement (\code{FALSE}) when using method `raref`.
#' @param multi Rarefy the table for multiple times, FALSE by default, indicate
#' the times of rarefaction want to be repeated, only validate for rarefaction.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref", depth = 1000, replace = TRUE,
#' multi = 3)
#' @return x An object of the class mina with @norm added.
#' @rdname norm_tab-mina
#' @exportMethod norm_tab

setMethod("norm_tab", signature("mina", "ANY", "ANY", "ANY", "ANY"),
          function(x, method, depth = 1000, replace = TRUE, multi = FALSE) {
             stop("Must specify a `method`. See `? norm_tab_method_list`")
          }
)

###############################################################################


#' @rdname norm_tab-mina
#' @exportMethod norm_tab

setMethod("norm_tab", signature("mina", "character", "ANY", "ANY", "ANY"),
          function(x, method, depth = 1000, replace = TRUE, multi = FALSE) {
              x@norm <- norm_tab(x@tab, method,
                                 depth = depth, replace = replace, multi = multi)
              x@norm <- as.matrix(x@norm)
              return(x)
          }
)

###############################################################################

#' Function for normalization, by total number of the reads in each sample.
#'
#' @param x A quantitative table with samples in columns and compositions in
#' rows.
#' @return A normalized quantitative table.
#' @keywords internal

norm_by_total <- function(x) {
    norm <- apply(x, 2, function(x) x / sum(x))
    return(norm)
}

###############################################################################

#' Function for normalization by rarefying the samples into the same depth,
#' modified from \pkg{phyloseq}, find it
#' \href{https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html}{here}.
#'
#' @param x A quantitative table with sample in columns and compositions in rows.
#' @param depth The depth for rarefying, 1000 by default.
#' @param replace Whether to sample with replacement (\code{TRUE}) or without
#' replacement (\code{FALSE}). Default \code{TRUE} for computational efficiency.
#' @return A normalized quantitative table.
#' @keywords internal

norm_by_raref <- function(x, depth = 1000, replace = TRUE) {
    # Make sure depth is of length 1.
    if (length(depth) > 1) {
        warning("`depth` had more than one value. ",
                "Using only the first. \n ... \n")
        depth <- depth[1]
    }

    if (depth <= 0) {
        stop("`depth` less than or equal to zero. ",
             "Need positive depth to work.")
    }

    # Remove samples contain few reads than `depth`
    if (min(colSums(x)) < depth) {
        rmsamples <- colnames(x)[colSums(x) < depth]
        message(length(rmsamples), " samples removed for low depth")
        x <- x[, colSums(x) >= depth]
    }

    if (ncol(x) <= 0) stop("No sample has more reads than `depth`.")

    x_norm <- apply(x, 2, rarefaction_subsample,
                    depth = depth, replace = replace)

    rownames(x_norm) <- rownames(x)
    #x_norm <- x_norm[rowSums(x_norm) > 0, ]
    return(x_norm)
}

###############################################################################

#' Rarefaction subsample function, one sample, modified from a internal function
#' in \pkg{phyloseq}, find it
#' \href{https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html}{here}.
#'
#' @importFrom methods as
#' @param x A column of quantitative table.
#' @param depth The depth for rarefying, 1000 by default.
#' @param replace Whether to sample with or without replacement, \code{TRUE} by
#' default for computational efficiency.
#' @keywords internal

rarefaction_subsample <- function(x, depth = 1000, replace = TRUE){
    # initial rarefied vector
    rare <- numeric(length(x))
    if (replace) {
        suppressWarnings(subsample <- sample(1:length(x), depth, replace = TRUE,
                                             prob = x))
    } else {
        df <- data.frame(sample = (1:length(x)), time = as.numeric(x))
        obs <- apply(df, 1, function(x) {
                         rep_len(x["sample"], x["time"])
        })
        obs <- unlist(obs, use.names = FALSE)
        suppressWarnings(subsample <- sample(obs, depth, replace = FALSE))
    }

    subsample_tab <- table(subsample)
    rare[as(names(subsample_tab), "integer")] <- subsample_tab
    return(rare)
}

###############################################################################

#' List of normalization methods supported in \code{\link[mina]{norm_tab}}
#'
#' Normalization methods should be specified by exact string match.
#' @format A list of character vectors.
#' \describe{
#'    \item{\code{raref}}{ By downsampling all samples to specific depth. }
#'    \item{\code{total}}{ Devided by the total read of each sample. }
#' }
#' @seealso \code{\link[mina]{norm_tab}}
#' @export
#' @examples
#' ? norm_tab_method_list

norm_tab_method_list <- list(
    raref = "raref",
    total = "total"
)
