#' Check the Tab and desTab file, used in setClass. Return TRUE if they are NULL or both quantitative and 
#' and descriptive files of same samples are included (i.e. the object is valid).
#'
#' @param object a mina object
#' @return TRUE if the object is valid or empty.

check_mina <- function(object){
    if (object@Tab == NULL & object@desTab == NULL) {
        return(TRUE)
        exit
    }
    errors <- character()
    errors <- c(errors, check_mina_qu(object))
    errors <- c(errors, check_mina_de(object))

    if (length(errors) == 0 ) TRUE else errors
}

#' Check the object and return TRUE if the object includes quantitative table.
#'
#' @param object a mina object.
#' @return TRUE if the object contains quantitative table and is not empty.
#' @examples
#' check_mina_qu(mina)
#' @export

check_mina_qu <- function(object){
    errors <- character()
    d <- dim(object@Tab)
    if (d[1] * d[2] == 0){
        msg <- paste0("Quantitative table is ", d[1], " * ", d[2], "!")
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}


#' Check the object and return TRUE if the object includes descriptive table.
#'
#' @param object a mina object.
#' @return TRUE if the object contains non-empty descriptive table and has the same samples as quantitative table.
#' @examples
#' check_mina_de(mina)
#' @export

check_mina_de <- function(object){
    errors <- character()
    nr <- nrow(object@desTab)
    nc <- ncol(object@desTab)

    if (nr * nc == 0){
        msg <- paste0("Descriptive table is ", nr, " * ", nc, "!")
        errors <- c(errors, msg)
    }

    samples1 <- sort(colnames(object@Tab))
    samples2 <- sort(object@desTab$SampleID)
    if (samples1 != samples2){
        msg <- "The samples in Tab and desTab are different!"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}
