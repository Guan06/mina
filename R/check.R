#' Check the object and return TRUE if the object includes quantitative table.
#'
#' @param mina a mina object.
#' @return TRUE if the object contains quantitative table.
#' @examples
#' check_mina_qu(mina)
#' export

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
#' @param mina a mina object.
#' @return TRUE if the object contains descriptive table.
#' @examples
#' check_mina_de(mina)
#' export

check_mina_de <- function(object){
    errors <- character()
    nr <- nrow(object@desTab)
    nc <- ncol(object@desTab)

    if (nr * nc == 0){
        msg <- paste0("Descriptive table is ", nr, " * ", nc, "!")
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}
