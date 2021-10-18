################################################################################

#' Calculate the network distance of `multi` and test the significance when
#' `perm` is defined.
#'
#' @importFrom stats dist
#' @param x An object of class `mina` with `multi` (and `perm` if sig is TRUE)
#' defined.
#' @param method The distance to be calculated, "spectra" and "Jaccard" are
#' available.
#' @param evk The first `evk` eigenvalues will be used for `spectra` distance,
#' the default is 100.
#' @param egv Wheather to output the eigenvectors for Spectral distance, the
#' defult is TRUE, only validate when `method == "spectra"`.
#' @param dir The folder to output the eigenvectors, only validate when `egv ==
#' TRUE`.
#' @param sig Whether to test the significance, if TRUE (by default), `perm` is
#' needed.
#' @param skip Whether to skip the comparison when the dimenstion of adjacency
#' matrix is smaller than setted `evk`.
#' @param ... Additional parameters.
#' @return x The same `mina` object with `net_dis` defined.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment")
#' maize <- net_dis(maize, method = "Jaccard")
#' @rdname net_dis-mina
#' @exportMethod net_dis

setMethod("net_dis", signature("mina", "ANY"),
          function(x, method, evk = 100, egv = TRUE, dir = "./",
                   sig = TRUE, skip = TRUE, ...) {
              stop("Must specify a `method`, see `? net_dis_method_list`.")
          }
)

###############################################################################

#' @rdname net_dis-mina
#' @exportMethod net_dis

setMethod("net_dis", signature("mina", "character"),
          function(x, method, evk = 100, egv = TRUE, dir = "./",
                   sig = TRUE, skip = TRUE, ...) {
              stopifnot(
                        method %in% c("spectra", "Jaccard"),
                        is.numeric(evk),
                        is.logical(c(egv, sig, skip)),
                        is.character(dir)
              )
              y_bs <- multi(x)
              y_pm <- perm(x)
              len <- length(y_bs)
              plen <- length(y_pm)

              if (sig) {
                  if (plen == 0) {
                      stop("Must have @perm for significance test!")
                  }else if (plen != len) {
                      stop("Must have same comparison numbers!")
                  }

                  dis_pm <- c()
              }

              dis_bs <- c()

              for (i in seq(1, len, 2)) {
                  this_m <- y_bs[[i]]
                  this_n <- y_bs[[i + 1]]

                  group_m <- names(y_bs)[i]
                  group_n <- names(y_bs)[i + 1]

                  ## calculate bootstrap distance
                  bs_len <- length(this_m)

                  if (method == "spectra") {
                      spectra_m <- spectra_n <- matrix(nrow = bs_len,
                                                       ncol = evk)
                      flag <- 0
                      for (j in 1 : bs_len) {
                          adj_m <- unlist(this_m[[j]])
                          adj_n <- unlist(this_n[[j]])
                          adj_m[is.na(adj_m)] <- 0
                          adj_n[is.na(adj_n)] <- 0

                          #skip this comparison if the adj_m and adj_n is
                          #smaller than evk
                          if (nrow(adj_m) < evk ||nrow(adj_n) < evk) {
                              flag <- 1
                              break
                          }
                          spectra_m[j, ] <- get_spectra(adj_m, k = evk)
                          spectra_n[j, ] <- get_spectra(adj_n, k = evk)
                       }
                      if (flag) next
                      seqs <- seq(1 : bs_len)
                      rownames(spectra_m) <- paste0(group_m, "_bs", seqs)
                      rownames(spectra_n) <- paste0(group_n, "_bs", seqs)
                      spectra_mn <- rbind(spectra_m, spectra_n)

                      if (egv) {
                          saveRDS(spectra_mn,
                                  file = paste0(dir, "/spectra_bs_", group_m,
                                                "_vs_", group_n, ".rds"))
                      }

                      this_dis_bs <- get_dis_df(dist(spectra_mn))

                      # filter intra group network comparison when comparing
                      # networks from different environments
                      if (group_m != group_n) {
                          r <- this_dis_bs$Group1 != this_dis_bs$Group2
                          this_dis_bs <- this_dis_bs[r, ]
                      }
                      dis_bs <- rbind(dis_bs, this_dis_bs)

                  } else if (method == "Jaccard") {
                      jaccard_mn <- c()
                      for (j1 in 1 : bs_len) {
                          adj_m <- unlist(this_m[[j1]])
                          adj_m[is.na(adj_m)] <- 0

                          m_j1 <- paste0(group_m, "_bs", j1)

                          for (j2 in 1 : bs_len) {
                              adj_n <- unlist(this_n[[j2]])
                              adj_n[is.na(adj_n)] <- 0

                              n_j2 <- paste0(group_n, "_bs", j2)

                              contrast <- sum(abs(adj_m - adj_n))
                              max <- sum(pmax(abs(adj_m), abs(adj_n)))

                              dis <- contrast / max

                              this <- data.frame(C1 = m_j1,
                                                 C2 = n_j2,
                                                 Distance = dis,
                                                 Group1 = group_m,
                                                 Group2 = group_n)
                              jaccard_mn <- rbind(jaccard_mn, this)
                          }
                      }
                      dis_bs <- rbind(dis_bs, jaccard_mn)
                  }

                  ## calculate permutation distance if sig == TRUE
                  if (sig) {
                      this_mp <- y_pm[[i]]
                      this_np <- y_pm[[i + 1]]

                      pm_len <- length(this_mp)

                      if (method == "spectra") {
                          spectra_mp <- spectra_np <- matrix(nrow = pm_len,
                                                             ncol = evk)

                          flag <- 0
                          for (k in 1 : pm_len) {
                              adj_mp <- unlist(this_mp[[k]])
                              adj_np <- unlist(this_np[[k]])
                              adj_mp[is.na(adj_mp)] <- 0
                              adj_np[is.na(adj_np)] <- 0

                              if (nrow(adj_mp) < evk ||nrow(adj_np) < evk) {
                                  flag <- 1
                                  break
                              }

                              spectra_mp[k, ] <- get_spectra(adj_mp, k = evk)
                              spectra_np[k, ] <- get_spectra(adj_np, k = evk)
                          }

                          if (flag) next
                          seqs <- seq(1 : pm_len)
                          rownames(spectra_mp) <- paste0(group_m, "_pm", seqs)
                          rownames(spectra_np) <- paste0(group_n, "_pm", seqs)
                          spectra_mnp <- rbind(spectra_mp, spectra_np)
                          if (egv) {
                              saveRDS(spectra_mnp,
                                      file = paste0(dir, "/spectra_pm_", group_m, "_vs_",
                                              group_n, ".rds"))
                          }
                          this_dis_pm <- get_dis_df(dist(spectra_mnp))

                          # filter intra group network comparison when
                          # comparing networks from different environments
                          if (group_m != group_n) {
                              r <- this_dis_pm$Group1 != this_dis_pm$Group2
                              this_dis_pm <- this_dis_pm[r, ]
                          }
                          dis_pm <- rbind(dis_pm, this_dis_pm)

                      } else if (method == "Jaccard") {
                          jaccard_mnp <- c()
                          for (k1 in 1 : pm_len) {
                              adj_mp <- unlist(this_mp[[k1]])
                              adj_mp[is.na(adj_mp)] <- 0

                              mp_k1 <- paste0(group_m, "_pm", k1)

                              for (k2 in 1 : pm_len) {
                                  adj_np <- unlist(this_np[[k2]])
                                  adj_np[is.na(adj_np)] <- 0

                                  np_k2 <- paste0(group_n, "_pm", k2)

                                  contrast <- sum(abs(adj_mp - adj_np))
                                  max <- sum(pmax(abs(adj_mp), abs(adj_np)))

                                  dis <- contrast / max
                                  this <- data.frame(C1 = mp_k1,
                                                     C2 = np_k2,
                                                     Distance = dis,
                                                     Group1 = group_m,
                                                     Group2 = group_n)

                                  jaccard_mnp <- rbind(jaccard_mnp, this)
                              }
                          }
                          dis_pm <- rbind(dis_pm, jaccard_mnp)
                      }
                  }
              }

              dis_bs(x) <- dis_bs

              ## get the stat table
              if (sig) {
                  dis_pm(x) <- dis_pm
                  dis_stat <- get_stat(dis_bs, dis_pm)
              } else {
                  dis_stat <- get_stat(dis_bs)
              }

              dis_stat(x) <- dis_stat
              return(x)
          }
)


################################################################################

#' Function for calculation of eigenvalue of given matrix.
#'
#' @importFrom RSpectra eigs_sym
#' @param x The input matrix.
#' @param k Get the first k eigenvalues.
#' @return y The vector of the first k eigenvalues.
#' @keywords internal

get_spectra <- function(x,  k = 100){
    x <- as.matrix(x)
    x[is.na(x)] <- 0
    spectra <- eigs_sym(x, k, opts = list(retvec = FALSE))
    y <- spectra$values
}

################################################################################

#' Function for getting distance data frame from `dist`.
#'
#' @importFrom reshape2 melt
#' @importFrom stringr str_detect
#' @param x The object of class `dist`.
#' @return The data frame of distance matrix.
#' @keywords internal

get_dis_df <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)] <- NA
    diag(x) <- NA
    x <- melt(x)

    colnames(x) <- c("C1", "C2", "Distance")
    x <- x[!is.na(x$Distance), ]

    if (str_detect(x$C1, "_bs")) {
        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C1), "_bs")))
        x <- x[, 1 : 4]

        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C2), "_bs")))
        x <- x[, 1 : 5]
    } else if (str_detect(x$C1, "_pm")) {
        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C1), "_pm")))
        x <- x[, 1 : 4]

        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C2), "_pm")))
        x <- x[, 1 : 5]
    }
    colnames(x)[4:5] <- c("Group1", "Group2")
    return(x)
}


###############################################################################

#' Function for distance statistic and significance test.
#'
#' @importFrom stats sd
#' @param x The bootstrap distance data frame.
#' @param p The permuation distance data frame.
#' @return The statistics of network comparison.
#' @keywords internal

get_stat <- function(x, p = NULL) {
    lst <- unique(x$Group1)
    len <- length(lst)

    grp_map <- c()
    y <- c()

    for ( i in 1 : len) {
        c1 <- lst[i]
        for (j in i : len) {
            c2 <- lst[j]
            d1 <- x[(x$Group1 == c1 & x$Group2 == c2), ]
            d2 <- x[(x$Group2 == c1 & x$Group1 == c2), ]
            ## skip the comparison skipped by spectra
            if (nrow(d1) == 0 & nrow(d2) == 0) next
            d <- rbind(d1, d2)
            d <- unique(d)

            if (i != j) d <- d[(d$Group1 != d$Group2), ]

            d$Compare <- paste0(c1, "_", c2)
            d <- d[, c("Compare", "Distance")]

            grp_map <- rbind(grp_map, c(paste0(c1, "_", c2), c1, c2))

            this_y <- data.frame(Compare = unique(d$Compare),
                                 Distance_Mean = mean(d$Distance),
                                 Distance_SD = sd(d$Distance))

            if (!is.null(p)) {
                dp1 <- p[(p$Group1 == c1 & p$Group2 == c2), ]
                dp2 <- p[(p$Group2 == c1 & p$Group1 == c2), ]
                #skip the comparisons skipped by spectra
                if (nrow(dp1) == 0 & nrow(dp2) == 0) next
                dp <- rbind(dp1, dp2)
                dp <- unique(dp)

                if (i != j) dp <- dp[(dp$Group1 != dp$Group2), ]
                dp$Compare <- paste0(c1, "_", c2)
                dp <- dp[, c("Compare", "Distance")]
                colnames(dp)[2] <- "Distance_PM"

                this_y$Distance_PM_Mean <- mean(dp$Distance_PM)
                this_y$Distance_PM_SD <- sd(dp$Distance_PM)

                d <- merge(d, dp)
                N <- nrow(d)
                this_y$N <- N
                this_y$p <- (sum(d$Distance_PM > d$Distance) + 1) / (N + 1)
            }

            y <- rbind(y, this_y)
        }
    }
    colnames(d) <- c("Compare", "Group1", "Group2")
    y <- merge(y, d)
    return(y)
}
