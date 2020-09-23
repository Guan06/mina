###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, samples present in both tables are remained. If @norm table exist in
#' the `mina` object, descriptive table will be filtered again to only keep
#' samples present in @norm.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class mina.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- fit_tabs(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' }
#' @return x The same object as input with fitted @tab, @des and @norm (if
#' defined).
#' @rdname fit_tabs-mina
#' @exportMethod fit_tabs

setMethod("fit_tabs", signature("mina"),
          function(x) {
              if (class(x@tab)[1] == "NULL" || class(x@des)[1] == "NULL") {
                  stop("Either @tab or @des of this object is missing!")
              }

              x@tab <- x@tab[rowSums(x@tab) > 0, ]

              samples1 <- as.character(colnames(x@tab))
              samples2 <- as.character(x@des$Sample_ID)

              x@des$Sample_ID <- as.factor(x@des$Sample_ID)

              inter <- intersect(samples1, samples2)
              x@tab <- x@tab[, colnames(x@tab) %in% inter]
              x@des <- x@des[x@des$Sample_ID %in% inter, ]

              # make the quantitative and descriptive files in the same order
              x@tab <- x@tab[rowSums(x@tab) > 0,
                             match(x@des$Sample_ID, colnames(x@tab))]

              ## filter the descriptive and quantitative again if @norm exists
              if (class(x@norm)[1] != "NULL") {
                  samples3 <- as.character(colnames(x@norm))
                  samples4 <- as.character(x@des$Sample_ID)
                  inter2 <- intersect(samples3, samples4)
                  x@norm <- x@norm[, colnames(x@norm) %in% inter2]
                  #x@tab <- x@tab[, colnames(x@tab) %in% inter2]
                  x@des <- x@des[x@des$Sample_ID %in% inter2, ]

                  # re-ordering
                  x@norm <- x@norm[rowSums(x@norm) > 0,
                                   match(x@des$Sample_ID, colnames(x@norm))]
                  #x@tab <- x@tab[rowSums(x@tab) > 0,
                  #               match(x@des$Sample_ID, colnames(x@tab))]
              }
              return(x)
          }
)
