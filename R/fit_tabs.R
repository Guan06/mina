###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, samples present in both tables are remained. If `norm` table exist in
#' the `mina` object, descriptive table will be filtered again to only keep
#' samples present in `norm`.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class mina.
#' @examples
#' {
#' data(maize)
#' maize <- fit_tabs(maize)
#' maize <- norm_tab(maize, method = "total")
#' maize <- fit_tabs(maize)
#' }
#' @return x The same object as input with fitted `tab`, `des` and `norm` (if
#' defined).
#' @rdname fit_tabs-mina
#' @exportMethod fit_tabs

setMethod("fit_tabs", signature("mina"),
          function(x) {
              if (class(tab(x))[1] == "NULL" || class(des(x))[1] == "NULL") {
                  stop("Either @tab or @des of this object is missing!")
              }

              tab(x) <- tab(x)[rowSums(tab(x)) > 0, ]

              samples1 <- as.character(colnames(tab(x)))
              samples2 <- as.character(des(x)$Sample_ID)

              des(x)$Sample_ID <- as.factor(des(x)$Sample_ID)

              inter <- intersect(samples1, samples2)
              tab(x) <- tab(x)[, colnames(tab(x)) %in% inter]
              des(x) <- des(x)[des(x)$Sample_ID %in% inter, ]

              # make the quantitative and descriptive files in the same order
              tab(x) <- tab(x)[rowSums(tab(x)) > 0,
                             match(des(x)$Sample_ID, colnames(tab(x)))]

              ## filter the descriptive and quantitative again if @norm exists
              if (class(norm(x))[1] != "NULL") {
                  samples3 <- as.character(colnames(norm(x)))
                  samples4 <- as.character(des(x)$Sample_ID)
                  inter2 <- intersect(samples3, samples4)
                  norm(x) <- norm(x)[, colnames(norm(x)) %in% inter2]
                  des(x) <- des(x)[des(x)$Sample_ID %in% inter2, ]

                  # re-ordering
                  norm(x) <- norm(x)[rowSums(norm(x)) > 0,
                                   match(des(x)$Sample_ID, colnames(norm(x)))]
              }
              return(x)
          }
)
