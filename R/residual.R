################################################################################

#' Function for unexplained variance ratio calculation indicated in
#' Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of
#' variance. Austral Ecology, 26: 32--46.
#'
#' @include all_classes.R all_generics.R
#' @param x An mina object with @dis and @des defined.
#' @param group The name(s) of column(s) defined as experimental setup group(s).
#'
#' @examples
#' x <- com_r2(x, group = "Soil")
#'
#' m <- new("mina", dis = distance_matrix, des = design_data_frame)
#' m <- com_r2(m)
#'
#' @return r2 The variance ratio cannot be explained by given groups.
#' @exportMethod com_r2

setMethod("com_r2", signature("mina", "ANY"), function(x, group) {
    stop("Must specify group(s) for unexplained variance ratio calculation!")
})

setMethod("com_r2", signature("mina", "character"), function(x, group) {
    r2 <- get_r2(x@dis, x@des, group = group)
    return(r2)
})


#' Function for unexplained variance ratio calculation indicated in
#' Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of
#' variance. Austral Ecology, 26: 32--46.
#'
#' @include all_classes.R all_generics.R
#'
#' @importFrom gdata drop.levels
#' @importFrom utils combn
#' @param x Dissimilarity / distance matrix which indicate variances.
#' @param des The descriptive table of samples which define the groups.
#' @param group The name(s) of column(s) used  as experimental setup group(s) in
#' descriptive file.
#'
#' @examples
#' r2 <- get_r2(x, des, group = c("Compartment", "Soil")
#'
#' @return r2 The variance ratio cannot be explained by given groups.
#' @export

get_r2 <- function(x, des, group = c("Host_genotype", "Compartment", "Soil",
                                     "Management")) {
    ## reformat the distance matrix x
    dis <- data.frame(t(combn(rownames(x), 2)), x[lower.tri(x)])
    colnames(dis) <- c("sample1", "sample2", "distance")
    dis <- dis[dis$distance != 0, ]

    len <- length(group)

    for (i in 1:len) {
        des$Group <- paste(des$Group, des[[group[i]]], sep="_")
    }

    des1 <- des[match(dis$sample1, des$Sample_ID), ]
    dis <- cbind(dis, des1$Group)

    des2 <- des[match(dis$sample2, des$Sample_ID), ]
    dis <- cbind(dis, des2$Group)

    colnames(dis)[4:5] <- c("sample1_group", "sample2_group")

    dis$distance_sq <- dis$distance ** 2
    SSt <- sum(dis$distance_sq) / nrow(x)

    dis_intra <- dis[dis$sample1_group == dis$sample2_group, ]
    lst <- levels(drop.levels(dis_intra$sample1_group))

    SSw <- 0
    for (l in lst) {
        this_group <- dis_intra[dis_intra$sample1_group == l, ]
        this_group <- drop.levels(this_group)
        n <- length(unique(c(levels(this_group$sample1), levels(this_group$sample2))))
        this_SSw <- sum(this_group$distance_sq) / n
        SSw <- SSw + this_SSw
    }

    r <- round(SSw / SSt, 3)
    return(r)
}
