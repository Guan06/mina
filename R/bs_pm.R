################################################################################

#' Inferring the network of different group of samples and test significance by
#' permutation.
#'
#' @include all_classes.R all_generics.R
#' @param x An object of class `mina` with @norm and @des defined.
#' @param group The column name of descriptive file @des for comparison.
#' @param g_size The cutoff of group size used for filtering, default is 88.
#' @param s_size The number of samples used for network inference during
#' bootstrap and permutation (when `sig` is TRUE), it should be smaller than
#' g_size / 2 to make sure the randomness; default is 30.
#' @param rm Filtering the components present in less than `per` of the samples
#' from compared groups, default TRUE.
#' @param per The percentage of present samples for filtering, default is 0.1.
#' @param sig Whether to test the significance, skip the permutation when set as
#' FALSE, default is TRUE.
#' @param bs The times for bootstrap network inference, default is 6.
#' @param pm The times for permuatated samples network inference, default is 6.
#' @param individual Whether to output the bootstrap and permutation results of
#' each comparison individually, default is FALSE.
#' @param out_dir The output directory if `individual` is TRUE, default is the
#' current working directory
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment", per = 0.5)
#' @return x The same object with @multi and @perm defined.
#' @rdname bs_pm-mina
#' @exportMethod bs_pm

setMethod("bs_pm", signature("mina", "ANY", "ANY", "ANY", "ANY", "ANY",
                             "ANY", "ANY", "ANY", "ANY", "ANY"),
          function(x, group, g_size = 88, s_size = 30, rm = TRUE, per = 0.1,
                   sig = TRUE, bs = 6, pm = 6,
                   individual = FALSE, out_dir = "./") {
              stop("Please specify a column in descriptive file for grouping
                   samples!")
          }
)

################################################################################

#' @rdname bs_pm-mina
#' @exportMethod bs_pm

setMethod("bs_pm", signature("mina", "character", "ANY", "ANY", "ANY", "ANY",
                            "ANY", "ANY", "ANY", "ANY", "ANY"),
          function(x, group, g_size = 88, s_size = 30, rm = TRUE, per = 0.1,
                   sig = TRUE, bs = 6, pm = 6,
                   individual = FALSE, out_dir = "./") {

              if (s_size >= g_size) {
                  stop("`s_size` can not be larger than `g_size`!")
              }

              mat <- x@norm
              des <- x@des

              mat <- mat[rowSums(mat) > 0, ]
              message(nrow(mat),
                      " components are used for bs_pm before filtering.")

              # fit the quantitative table with descriptive table
              fit <- intersect(colnames(mat), des$Sample_ID)
              mat <- mat[, colnames(mat) %in% fit]
              des <- des[des$Sample_ID %in% fit, ]
              message(nrow(des), " samples are used for bs_pm before filtering.")

              # filter the group does not have enough samples (i.e. < g_size)
              des[[group]] <- as.factor(des[[group]])
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
              if (!individual){
                  index <- 1
                  y_bs <- list()
                  y_pm <- list()
              }

              for (m in 1 : len) {
                  group_m <- lst[m]
                  des_m <- des[des[[group]] == group_m, ]
                  mat_m <- mat[, colnames(mat) %in% des_m$Sample_ID]
                  num_m <- nrow(des_m)

                  for (n in m : len) {
                      if (n == m) {
                          group_n <- group_m
                          this_mat_m <- this_mat_n <- mat_m
                          mat_mn <- mat_n <- mat_m
                          num_mn <- num_n <- num_m

                          # filter if rm is TRUE
                          if (rm) {
                              mat_mn <- filter_mat(mat_mn, p = s_size * per)
                              this_mat_m <- this_mat_n <- mat_mn
                              num_mn <- ncol(mat_mn)
                              if (num_mn < s_size) {
                                  stop("Not enough samples for ", group_m,
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
                              mat_mn <- filter_mat(mat_mn, p = size * per * 2)
                              this_mat_m <- mat_mn[, colnames(mat_mn) %in%
                                              des_m$Sample_ID]
                              this_mat_n <- mat_mn[, colnames(mat_mn) %in%
                                              des_n$Sample_ID]

                              num_m <- ncol(this_mat_m)
                              num_n <- ncol(this_mat_n)
                              num_mn <- num_m + num_n

                              if (num_m < s_size || num_n < s_size) {
                                  stop("Not enough samples for ", group_m,
                                       " v.s. ", group_n,
                                       " bootstrap after filtering!")
                              }
                          } else {
                              this_mat_m <- mat_m
                              this_mat_n <- mat_n
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

                          mat_bs_m <- this_mat_m[, bs_m]
                          mat_bs_n <- this_mat_n[, bs_n]

                          cor_m <- adj(mat_bs_m, method = "spearman")
                          cor_n <- adj(mat_bs_n, method = "spearman")

                          cor_m[is.na(cor_m)] <- 0
                          cor_n[is.na(cor_n)] <- 0

                          MLST[[b]] <- cor_m
                          NLST[[b]] <- cor_n
                          names(MLST)[b] <- paste0(group_m, "_", b)
                          names(NLST)[b] <- paste0(group_n, "_", b)
                      }

                      if (individual) {
                          prefix <- paste0(out_dir, group_m, "_vs_", group_n)
                          saveRDS(list(MLST), paste0(prefix, "_bs1.rds"))
                          saveRDS(list(NLST), paste0(prefix, "_bs2.rds"))
                      } else{
                          y_bs[index] <- list(MLST)
                          y_bs[(index + 1)] <- list(NLST)
                          names(y_bs)[index : (index + 1)] <- c(group_m, group_n)
                      }

                      rm(MLST, NLST)
                      gc(reset = TRUE)

                      # start permutation
                      if (sig) {
                          MPLST <- list()
                          NPLST <- list()

                          for (p in 1 : pm) {
                              pm_mn <- sample.int(num_mn, s_size * 2)
                              pm_m <- pm_mn[1 : s_size]
                              pm_n <- pm_mn[(s_size + 1) : (s_size * 2)]

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

                          if (individual) {
                              prefix <- paste0(out_dir, group_m, "_vs_",
                                               group_n)
                              saveRDS(list(MPLST), paste0(prefix, "_pm1.rds"))
                              saveRDS(list(NPLST), paste0(prefix, "_pm2.rds"))
                          } else {
                              y_pm[index] <- list(MPLST)
                              y_pm[(index + 1)] <- list(NPLST)

                              names(y_pm)[index] <- group_m
                              names(y_pm)[index + 1] <- group_n
                          }
                          rm(MPLST, NPLST)
                          gc(reset = TRUE)
                      }
                      if (!individual) index <- index + 2
                  }
              }
              if(!individual) {
                  x@multi <- y_bs
                  if (sig) x@perm <- y_pm
              }
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
#' @keywords internal

filter_mat <- function(x, p) {
    x_bi <- x
    x_bi[x_bi > 0] <- 1
    x <- x[rowSums(x_bi) > p, ]
    x <- x[, colSums(x) > 0]
    return(x)
}
