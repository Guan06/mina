###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, samples present in both tables are remained. If @norm table exist in
#' the `mina` object, descriptive table will be filtered again to only keep
#' samples present in @norm.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class mina.
#' @examples
#' data(maize)
#' maize <- fit_tabs(maize)
#' @return x The same object as input with fitted @tab and @des.
#' @exportMethod fit_tabs

setMethod("fit_tabs", signature("mina"),
          function(x) {
              if (class(x@tab) == "NULL" || class(x@des) == "NULL") {
                  stop("Either @tab or @des of this object is missing!")
              }

              samples1 <- as.character(colnames(x@tab))
              samples2 <- as.character(x@des$Sample_ID)

              x@des$Sample_ID <- as.factor(x@des$Sample_ID)

              inter <- intersect(samples1, samples2)
              x@tab <- x@tab[, colnames(x@tab) %in% inter]
              x@des <- x@des[x@des$Sample_ID %in% inter, ]
              #x@des <- x@des %>% filter(Sample_ID %in% inter)

              ## filter the descriptive again if @norm exists
              if (class(x@norm) != "NULL") {
                  samples3 <- as.character(colnames(x@norm))
                  samples4 <- as.character(colnames(x@des$Sample_ID))
                  inter2 <- intersect(samples3, samples4)
                  x@norm <- x@norm[, colnames(x@norm) %in% inter2]
                  x@des <- x@des[x@des$Sample_ID %in% inter2, ]
                  #x@des <- x@des %>% filter(Sample_ID %in% inter2)
              }
              return(x)
          }
)
