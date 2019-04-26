###############################################################################

#' Filter the quantitative and descriptive table to make them have the same
#' samples, samples present in both tables are remained.
#'
#' @include all_classes.R all_generics.R
#' @import dplyr
#' @param x An object of class mina.
#' @examples
#' x <- fit_tabs(x)
#' @return x The same object as input with fitted @tab and @des_tab.
#' @exportMethod fit_tabs

setMethod("fit_tabs", signature("ANY"),
          function(x) {
              stop("You must specify a `mina` class object as input.")
          }
)

setMethod("fit_tabs", signature("mina"),
          function(x) {
              if (class(x@tab) == "NULL" || class(x@des_tab) == "NULL") {
                  stop("Either @tab or @des_tab of this object is missing!")
              }

              samples1 <- as.character(sort(colnames(x@tab)))
              samples2 <- as.character(sort(x@des_tab$Sample_ID))

              inter <- intersect(samples1, samples2)
              x@tab <- x@tab[, colnames(x@tab) %in% inter]
              x@des_tab <- x@des_tab %>% filter(Sample_ID %in% inter)

              return(x)
          }
)
