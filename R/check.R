###############################################################################

#' Check the tab and des_tab file. Return TRUE if they are NULL or both
#' quantitative and descriptive files of same samples are included (i.e. the
#' object is valid).
#'
#' @param x An object of class mina.
#' @return TRUE if the object is valid.
#' @examples
#' check_mina(x)
#' @export

check_mina <- function(x){
    if (class(x@tab) == "NULL" && class(x@des_tab) == "NULL") {
        stop("An empty (neither @tab or @des_tab) object of the class mina!")
    }

    errors <- character()
    if (!check_mina_qu(x)) errors <- c(errors, check_mina_qu(x))
    if (!check_mina_de(x)) errors <- c(errors, check_mina_de(x))

    if (length(errors) == 0 ) TRUE else stop(errors)
}

###############################################################################

#' Check the object and return TRUE if the object includes quantitative table.
#'
#' @param x An object of class mina.
#' @return TRUE if the object contains quantitative table and is not empty.
#' @examples
#' check_mina_qu(x)
#' @export

check_mina_qu <- function(x){
    errors <- character()
    if (class(x@tab) == "NULL") stop("The @tab of this object does not exist!")
    #d <- dim(x@tab)
    #message (paste0("The @tab is ", d[1], " * ", d[2], "."))
    TRUE
}

###############################################################################

#' Check the object and return TRUE if the object includes descriptive table
#' contains the same samples as quantitative table.
#'
#' @param x An object of class mina.
#' @return TRUE if the object contains non-empty descriptive table and has the
#' same samples as quantitative table.
#' @examples
#' check_mina_de(x)
#' @export

check_mina_de <- function(x){
    errors <- character()

    if (!check_mina_qu(x)) errors <- c(errors, check_mina_qu(x))

    if(class(x@des_tab) == "NULL") {
        stop("The @des_tab of this object does not exist!")
    }

    samples1 <- as.character(sort(colnames(x@tab)))
    samples2 <- as.character(sort(x@des_tab$Sample_ID))

    if (!identical(samples1, samples2)){
        msg <- "The samples in @tab and @des_tab are different!"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else stop(errors)
}