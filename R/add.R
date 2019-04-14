#' Add tables to the object
#'
#' This includes functions for adding quantitative table (e.g. OTU table)
#' to mina object for later analysis.
#'
#' @param mina An object of the class mina.
#' @param file The quantitative table.
#' @return mina with quantitative table.
#'
#' @examples
#' addTab(mina, "otu_table.txt")
#'
#' @export

addTab <- function(mina, file){
    tab <- read.table(file, header=T, sep="\t", check.names=F)


    return(mina)
}

#' Add tables to the object
#'
#' This includes functions for adding descriptive table (e.g. design file)
#' to mina object for later analysis.
#'
#' @param mina An object of the class mina.
#' @param file The descriptive table.
#' @return mina with descriptive table.
#'
#' @examples
#' addTab(mina, "design.txt")
#'
#' @export

addDes <- function(mina, file){
    tab <- read.table(file, header=T, sep="\t")
    if (colnames(tab)[1] != "SampleID") {
        print("The first column of descriptive file should be SampleID!")
        exit
    }


    return(mina)
}
