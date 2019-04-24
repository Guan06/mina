################################################################################
#' Add quantitative table to the object.
#'
#' This includes functions for adding quantitative table (e.g. OTU table) to
#' mina object for later analysis.
#'
#' @param x An object of the class mina.
#' @param file The table to be added to the object.
#' @return x An object of the class mina with @tab added.
#'
#' @examples
#' add_tab(x, "otu_table.txt")
#' add_tab(x, "asv_table.rds")
#'
#' @export

add_tab <- function(x, file) {
    type <- unlist(strsplit(file, ".", fixed = TRUE))[2]

    if (type == "rds") {
        tab <- readRDS(file)
    }else if (type == "txt") {
        tab <- read.table(file, header=T, sep="\t", check.names=F)
    }else{
        stop("Unknown input file format!")
    }

    tab <- tab[rowSums(tab)>0, colSums(tab) >0]
    x@tab <- tab

    return(x)
}

################################################################################
#' Add descriptive table to the object.
#'
#' This includes functions for adding descriptive table (e.g. design file)
#' to mina object for later analysis. The samples in first column "Sample_ID"
#' should be the same as the samples in the quantitative table.
#'
#' @param x An object of the class mina with @tab defined.
#' @param file The descriptive table.
#' @return x An object of the class mina with @des_tab added.
#'
#' @examples
#' add_des(x, "design.txt")
#'
#' @export

add_des <- function(x, file) {
    des_tab <- read.table(file, header=T, sep="\t")

    if (colnames(tab)[1] != "Sample_ID") {
        stop("The first column of descriptive file should be Sample_ID!")
    }

    x@des_tab <- des_tab

    c <- check_mina_de(x)

    if (c == TRUE) x else c
}
