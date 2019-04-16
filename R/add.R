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
    tab <- tab[rowSums(tab)>0, colSums(tab) >0]
    mina@Tab <- tab

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
#' @return mina with descriptive table.
#'
#' @examples
#' addDes(mina, "design.txt")
#'
#' @export

addDes <- function(mina, file){
    des.tab <- read.table(file, header=T, sep="\t")
    if (colnames(tab)[1] != "SampleID") {
        print("The first column of descriptive file should be SampleID!")
        exit
    }

    mina@desTab <- des.tab

    c <- check_mina_de(mina)

    if (c == TRUE) mina else c
}
