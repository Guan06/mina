###############################################################################

#' Normalize the quantitative table with matrix input.
#'
#' @include all_classes.R all_generics.R
#' @param x A matrix of the quantitative table.
#' @param method The method used for normalization.
#' @param depth (optional) The depth for subsampling by rarefying, using the
#' minimum sample depth by default.
#' @param replace (optional) Whether to sample with replacement (\code{TRUE} by
#' default) or without replacement (\code{FALSE}) when using method `raref`.
#' @examples
#' norm_tab(x, method = "raref", depth = 1000, replace = TRUE)
#' norm_tab(x, method = "total")
#' @return norm_x Normalized matrix of the quantitative table.
#' @exportMethod norm_tab

setMethod("norm_tab", signature("matrix", "character", "ANY", "ANY"),
          function(x, method, depth = 1000, replace = TRUE) {
              if (method == "raref") {
                  norm_x <- norm_by_raref(x, depth = depth, replace = replace)
              }
              if (method == "total") {
                  norm_x <- norm_by_total(x)
              }
              return(norm_x)
          }
)

#' Normalize the quantitative table with mina input.
#'
#' @include all_classes.R all_generics.R
#'
#' @param x An object of the class mina with @tab defined.
#' @param method The method used for normalization.
#' @param depth (optional) The depth for subsampling by rarefying, 1000 by
#' default.
#' @param replace (optional) Whether to sample with replacement (\code{TRUE} by
#' default) or without replacement (\code{FALSE}) when using method `raref`.
#' @examples
#' x <- norm_tab(x, method = "raref", depth = 1000, replace = TRUE)
#' x <- norm_tab(x, method = "total")
#' @return x An object of the class mina with @norm added.
#' @exportMethod norm_tab

setMethod("norm_tab", signature("mina", "ANY", "ANY", "ANY"),
          function(x, method, depth = 1000, replace = TRUE) {
             stop("Must specify a `method`. See `? norm_tab_method_list`")
          }
)

setMethod("norm_tab", signature("mina", "character", "ANY", "ANY"),
          function(x, method, depth = 1000, replace = TRUE) {
              x@norm <- norm_tab(x@tab, method,
                                 depth = depth, replace = replace)
              return(x)
          }
)

###############################################################################

#' Function for normalization, by total number of the reads in each sample.
#'
#' @param x A quantitative table with samples in cloumns and compostions in
#' rows.
#' @examples
#' norm_tab <- norm_by_total(mina@tab)
#' @return A normalized quantitative table.
#' @keywords internal

norm_by_total <- function(x) {
    norm <- apply(x, 2, function(x) x / sum(x))
    return(norm)
}

###############################################################################

#' Function for normalization by rarefying the samples into the same depth,
#' modified from \code{\link[phyloseq]{rarefy_even_depth}}.
#'
#' @param x A quantitative table with sample in columns and compostions in rows.
#' @param depth (optional) The depth for rarefying, 1000 by default.
#' @param replace (optional) Whether to sample with replacement (\code{TRUE}) or
#' without replacement (\code{FALSE}). Default \code{TRUE} for computational
#' efficiency.
#' @examples
#' norm_tab <- norm_by_raref(x, depth = 1000, replace = FALSE)
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
    x_norm <- x_norm[rowSums(x_norm) > 0, ]
    return(x_norm)
}

#' Rarefaction subsample function, one sample, modified from
#' \code{\link[phyloseq]{rarefaction_subsample}} which is a internal function.
#' Resample with replacement for default. As mentioned in
#' \pkg{phyloseq}, this is set for computational efficiency.
#'
#' @param x A column of quantitative table.
#' @param depth The depth for rarefying.
#' @param replace Whether to sample with or without replacement.
#' @examples
#' rarefaction_subsample(otu_table[, 1], depth = 1000, replace = FALSE)
#' @keywords internal

rarefaction_subsample <- function(x, depth, replace = TRUE){
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
#'    \item{raref By downsampling all samples to specific depth.}
#'    \item{total Devided by the total read of each sample.}
#' }
#' @seealso \code{\link[mina]{norm_tab}}
#' @export
#' @examples
#' ? norm_tab_method_list

norm_tab_method_list <- list(
    raref = "raref",
    total = "total"
)
