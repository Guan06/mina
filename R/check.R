###############################################################################

#' Check the tab and des file. Return TRUE if they are NULL or both
#' quantitative and descriptive files of same samples are included (i.e. the
#' object is valid).
#'
#' @param x An object of class mina.
#' @return TRUE if the object is valid.
#' @examples
#' data(maize)
#' check_mina(maize)
#' @export

check_mina <- function(x) {
    if (class(x@tab)[1] == "NULL" && class(x@des)[1] == "NULL") {
        stop("An empty (neither @tab or @des) object of the class mina!")
    }

    errors <- character()
    if (!check_mina_qu(x)) errors <- c(errors, check_mina_qu(x))
    if (!check_mina_de(x)) errors <- c(errors, check_mina_de(x))

    if (length(errors) == 0 ) TRUE else stop(errors)
}

###############################################################################

#' Check the object and return TRUE if the object includes quantitative table.
#'
#' @param x An object of class mina with @tab defined.
#' @return TRUE if the object contains quantitative table and is not empty.
#' @examples
#' data(maize)
#' check_mina_qu(maize)
#' @export

check_mina_qu <- function(x) {
    errors <- character()
    if (class(x@tab)[1] == "NULL") stop("The @tab of this object does not exist!")
    TRUE
}

###############################################################################

#' Check the object and return TRUE if the object includes descriptive table
#' contains the same samples as quantitative table.
#'
#' @param x An object of class mina with @tab and @des defined.
#' @return TRUE if the object contains non-empty descriptive table and has the
#' same samples as quantitative table.
#' @examples
#' data(maize)
#' check_mina_de(maize)
#' @export

check_mina_de <- function(x) {
    errors <- character()

    if (!check_mina_qu(x)) errors <- c(errors, check_mina_qu(x))

    if(class(x@des)[1] == "NULL") {
        stop("The @des of this object does not exist!")
    }

    samples1 <- as.character(sort(colnames(x@tab)))
    samples2 <- as.character(sort(x@des$Sample_ID))

    if (!identical(samples1, samples2)){
        msg <- "The samples in @tab and @des are different!"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else stop(errors)
}
