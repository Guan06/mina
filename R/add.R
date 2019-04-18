#' Add tables to the object
#'
#' This includes functions for adding quantitative table (e.g. OTU table)
#' to mina object for later analysis.
#'
#' @param mina An object of the class mina.
#' @param file The quantitative table.
#' @return mina An object of the class mina with quantitative table added.
#'
#' @examples
#' add_tab(mina, "otu_table.txt")
#'
#' @export

add_tab <- function(mina, file) {
    tab <- read.table(file, header=T, sep="\t", check.names=F)
    tab <- tab[rowSums(tab)>0, colSums(tab) >0]
    mina@tab <- tab

    return(mina)
}

#' Add tables to the object
#'
#' This includes functions for adding descriptive table (e.g. design file)
#' to mina object for later analysis. The samples in first column "SampleID" 
#' should be the same as the samples in the quantitative table.
#'
#' @param mina An object of the class mina with Tab defined.
#' @param file The descriptive table.
#' @return mina An object of the class mina with descriptive table added.
#'
#' @examples
#' add_des(mina, "design.txt")
#'
#' @export

add_des <- function(mina, file) {
    des.tab <- read.table(file, header=T, sep="\t")

    if (colnames(tab)[1] != "SampleID") {
        stop("The first column of descriptive file should be SampleID!")
    }

    mina@des_tab <- des.tab

    c <- check_mina_de(mina)

    if (c == TRUE) mina else c
}

#' Add normalized quantitative table to object mina
#'
#' The normalized quantitative table should contains the exactly same samples
#' as quantitative table and descriptive table in the same mina object.
#'
#' @param mina An object of the class mina.
#' @param norm The normalized quantitative table or the file contains the table.
#'
#' @examples
#' add_norm(mina, norm)
#' add_norm(mina, "otu_table_norm.txt")
#' @return mina An object of the class mina with normalized quantitative table
#' added.
#'
#' @export

add_norm <- function(mina, norm) {
    if (class(norm) == "character") {
        norm <- read.table(norm, header=T, sep="\t", check.names=F)
    }
    norm <- norm[rowSums(norm) > 0, colSums(norm) > 0]
    mina@norm <- norm
    return (mina)
}
