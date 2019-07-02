################################################################################

#' Inferring the network of different group of samples and test significance by
#' permutation.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class `mina` with @norm and @des_tab defined.
#' @param group The column name of descriptive file @des_tab for comparison.
#' @param g_size The cutoff of group size used for filtering, default 88.
#' @param s_size The number of samples used for network inference during
#' bootstrap and permutation (when sig == TRUE), it should be smaller than
#' g_size / 2 to make sure the randomness; default 30.
#' @param rm Filtering the components present in less than 20% of the samples
#' from compared groups, default TRUE.
#' @param sig Whether to test the significance, skip the permutation when sig ==
#' FALSE, default TRUE.
#' @param bs The times for bootstrap network inference, default 6
#' @param pm The times for permuatated samples network inference, default 6.
#' @examples
#' x <- bs_pm(x, group = "Compartment")
#' x <- bs_pm(x, group = "Compartment", g_size = 100, s_size = 50, rm = TRUE,
#' sig = TRUE, bs = 6, pm = 6)')
#' @return x The same object with @multi and @perm defined.
#' @exportMethod bs_pm

setMethod("bs_pm", signature("mina", "ANY", "numeric", "numeric", "logical",
                             "logical", "numeric", "numeric"),
          function(x, group, g_size = 88, s_size = 30, rm = TRUE, sig = TRUE,
                   bs = 6, pm = 6) {
              stop("Please specify a column in descriptive file for grouping
                   samples!")
          }
)

setMethod("bs_pm", signature("mina", "character", "numeric", "numeric", "logical",
                            "logical", "numeric", "numeric"),
          function(x, group, g_size = 100, s_size = 50, rm = TRUE, sig = TRUE,
                   bs = 6, pm = 6) {

              if (s_size >= g_size) {
                  stop("`s_size` can not be larger than `g_size`!")
              }

              mat <- x@norm
              des <- x@des_tab

              # fit the quantitative table with descriptive table
              fit <- intersect(colnames(mat), des$Sample_ID)
              mat <- mat[, colnames(mat) %in% fit]
              des <- des[des$Sample_ID %in% fit, ]
              message(nrow(des), " samples are used for bs_pm before filtering.")

              # filter the group does not have enough samples (i.e. < g_size)
              lst <- levels(des[[group]])
              len <- length(lst)

              for (i in 1 : len) {
                  g <- lst[i]
                  size <- sum(des[[group]] == g)
                  if (size < g_size) lst[i] <- NA
              }

              lst <- lst[!is.na(lst)]
              len <- length(lst)

              des <- des[des[[group]] %in% lst, ]
              mat <- mat[, colnames(mat) %in% des$Sample_ID]
              message(len, " groups with ", ncol(mat),
                      " samples used for bootstrap.")

              # start bootstrap and permutation
              index <- 1
              y_bs <- list()
              y_pm <- list()

              for (m in 1 : len) {
                  group_m <- lst[m]
                  des_m <- des[des[[group]] == group_m, ]
                  mat_m <- mat[, colnames(mat) %in% des_m$Sample_ID]
                  num_m <- nrow(des_m)

                  for (n in m : len) {
                      if (n == m) {
                          group_n <- group_m
                          mat_mn <- mat_n <- mat_m
                          num_mn <- num_n <- num_m

                          # filter if rm is TRUE
                          if (rm) {
                              mat_mn <- filter_mat(mat_mn, p = s_size * 0.1)
                              num_mn <- ncol(mat_mn)
                              if (num_mn < s_size) {
                                  stop("Not enough samples for", group_m,
                                       " bootstrap after filtering!")
                              }
                          }

                      } else {
                          group_n <- lst[n]
                          des_n <- des[des[[group]] == group_n, ]
                          mat_n <- mat[, colnames(mat) %in% des_n$Sample_ID]
                          num_n <- nrow(des_n)
                          mat_mn <- cbind(mat_m, mat_n)
                          num_mn <- num_m + num_n

                          # filter if rm is TRUE
                          if (rm) {
                              mat_mn <- filter_mat(mat_mn, p = size * 0.2)
                              mat_m <- mat_mn[, colnames(mat_mn) %in%
                                              des_m$Sample_ID]
                              mat_n <- mat_mn[, colnames(mat_mn) %in%
                                              des_n$Sample_ID]

                              num_m <- ncol(mat_m)
                              num_n <- ncol(mat_n)
                              num_mn <- num_m + num_n

                              if (num_m < s_size || num_n < s_size) {
                                  stop("Not enough samples for ", group_m,
                                       " v.s. ", group_n,
                                       " bootstrap after filtering!")
                              }
                          }
                      }

                      # re-normalization
                      mat_mn <- apply(mat_mn, 2, function(x) x / sum(x))

                      # start bootstrap
                      MLST <- list()
                      NLST <- list()

                      for (b in 1 : bs) {
                          bs_m <- sample.int(num_m, s_size)
                          bs_n <- sample.int(num_n, s_size)

                          mat_bs_m <- mat_m[, bs_m]
                          mat_bs_n <- mat_n[, bs_n]

                          cor_m <- adj(mat_bs_m, method = "spearman")
                          cor_n <- adj(mat_bs_n, method = "spearman")

                          cor_m[is.na(cor_m)] <- 0
                          cor_n[is.na(cor_n)] <- 0

                          MLST[[b]] <- cor_m
                          NLST[[b]] <- cor_n
                          names(MLST)[b] <- paste0(group_m, "_", b)
                          names(NLST)[b] <- paste0(group_n, "_", b)
                      }

                      #y_bs[index] <- list(MLST, NLST)
                      #names(y_bs)[index] <- paste0(group_m, "_", group_n)
                      y_bs[index] <- list(MLST)
                      y_bs[(index + 1)] <- list(NLST)
                      names(y_bs)[index : (index + 1)] <- c(group_m, group_n)

                      rm(MLST, NLST)
                      gc(reset = T)

                      # start permuationg
                      if (sig) {
                          MPLST <- list()
                          NPLST <- list()

                          for (p in 1 : pm) {
                              pm_mn <- sample.int(num_mn, s_size * 2)
                              pm_m <- pm_mn[1 : 50]
                              pm_n <- pm_mn[51 : 100]

                              mat_pm_m <- mat_mn[, pm_m]
                              mat_pm_n <- mat_mn[, pm_n]

                              cor_pm <- adj(mat_pm_m, method = "spearman")
                              cor_pn <- adj(mat_pm_n, method = "spearman")

                              cor_pm[is.na(cor_pm)] <- 0
                              cor_pn[is.na(cor_pn)] <- 0

                              MPLST[[p]] <- cor_pm
                              NPLST[[p]] <- cor_pn
                              names(MPLST)[p] <- paste0(group_m, "_", p)
                              names(NPLST)[p] <- paste0(group_n, "_", p)
                          }

                          #y_pm[index] <- list(MPLST, NPLST)
                          #names(y_pm)[index] <- paste0(group_m, "_", group_n)
                          y_pm[index] <- list(MPLST)
                          y_pm[(index + 1)] <- list(NPLST)

                          names(y_pm)[index] <- group_m
                          names(y_pm)[index + 1] <- group_n

                          rm(MPLST, NPLST)
                          gc(reset = T)
                      }
                      index <- index + 2
                  }
              }
              x@multi <- y_bs
              if (sig) x@perm <- y_pm
              return(x)
          }
)

###############################################################################

#' Function for filtering of matrix, rows present in less than `p` columns will
#' be removed. After row filtering, columns with a sum of 0 wil be removed.
#'
#' @param x The input matrix to be filtered.
#' @param p The cutoff for non-zero column number.
#' @return x The same matrix after filtering.
#' @examples
#' x <- filter_mat(x, p = 10)
#' @keywords internal

filter_mat <- function(x, p) {
    x_bi <- x
    x_bi[x_bi > 0] <- 1
    x <- x[rowSums(x_bi) > p, ]
    x <- x[, colSums(x) > 0]
    return(x)
}
