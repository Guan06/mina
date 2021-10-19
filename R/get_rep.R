###############################################################################

#' Get the representative community members.
#'
#' @include all_classes.R all_generics.R
#' @param x A quantitative matrix with samples in columns and compositions in
#' rows.
#' @param top The percent of the most abundant and prevalent members.
#' @return The matrix with samples in columns and representative compositions in
#' rows.
#' @examples
#' data(maize_asv)
#' maize_asv_rep <- get_rep(maize_asv, top = 5)
#' @rdname get_rep-matrix
#' @exportMethod get_rep

setMethod("get_rep", signature("matrix", "ANY"),
          function(x, top = 5) {
    stopifnot(is.numeric(top))
    asv <- x
    ## get the most abundant members
    ra <- rowSums(asv) / ncol(asv)
    ra <- ra[order(-ra)]
    top_num <- round(nrow(asv) * top / 100)
    ra_top <- names(ra[1 : top_num])

    ## get the most prevalent members
    occu <- rowSums(asv > 0) / ncol(asv)
    occu <- occu[order(-occu)]
    top_num <- round(nrow(asv) * top / 100)
    occu_top <- names(occu[1 : top_num])

    rep <- intersect(ra_top, occu_top)
    rep_asv <- asv[rownames(asv) %in% rep, ]
    print(dim(rep_asv))
    
    rep_asv <- apply(rep_asv, 2, function(x) x / sum(x))
    }
)
###############################################################################

#' Get the representative community members.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of the class `mina` with @norm define.
#' @param top The percent of the most abundant and prevalent members.
#' @return The same object with @norm replaced by the representative members.
#' @examples
#' maize <- new("mina", tab = maize_asv, des = maize_des)
#' maize <- get_rep(maize, top = 5)
#' @rdname get_rep-mima
#' @exportMethod get_rep

setMethod("get_rep", signature("mina", "ANY"),
          function(x, top = 5) {
    
              stopifnot(is.numeric(top))
              norm(x) <- get_rep(norm(x))
              return(x)
    }
)
