################################################################################

#' Calculate the network distance of bootstrap and permutation when appliable.
#'
#' @importFrom stats dist
#' @param x The folder store the network inference results.
#' defined.
#' @param method The distance to be calculated, "spectra" and "Jaccard" are
#' available.
#' @param evk The first `evk` eigenvalues will be used for `spectra` distance,
#' the default is 100.
#' @param sig Whether to test the significance, if TRUE (by default),
#' permutation results should be included in the folder `x`.
#' @param skip Whether to skip the comparison when the dimenstion of adjacency
#' matrix is smaller than setted `evk`.
#' @return x The `mina` object with @dis_bs, @dis_pm and @dis_stat.
#' @examples
#' \dontrun{
#' data(maize)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment", individual = TRUE, out_dir =
#' "./individual_bs_pm/")
#' maize_stat1 <- net_dis_indi("./individual_bs_pm/", method = "spectra")
#' maize_stat2 <- net_dis_indi("./individual_bs_pm/", method = "Jaccard")
#' maize_stat3 <- net_dis_indi("./individual_bs_pm/", method = "spectra",
#' evk = 100, skip = TRUE)
#' }
#' @rdname net_dis_indi
#' @exportMethod net_dis_indi

setMethod("net_dis_indi", signature("character", "ANY", "ANY", "ANY", "ANY"),
          function(x, method, evk = 100, sig = TRUE, skip = TRUE) {
              stop("Must specify a `method`, see `? net_dis_method_list`.")
          }
)

###############################################################################

#' @inheritParams net_dis_indi
#' @rdname net_dis_indi
#' @exportMethod net_dis_indi

setMethod("net_dis_indi", signature("character", "character",
                                    "ANY", "ANY", "ANY"),
          function(x, method, evk = 100, sig = TRUE, skip = TRUE) {
              bs1_files <- sort(list.files(x, pattern = "_bs1.rds",
                                           full.names = TRUE))
              bs2_files <- sort(list.files(x, pattern = "_bs2.rds",
                                           full.names = TRUE))
              dis_bs <- c()

              if (sig){
                  pm1_files <- sort(list.files(x, pattern = "_pm1.rds",
                                               full.names = TRUE))
                  pm2_files <- sort(list.files(x, pattern = "_pm2.rds",
                                               full.names = TRUE))
                  dis_pm <- c()
              }

              log <- c()

              len <- length(bs1_files)

              for (i in 1 : len) {
                  bs1 <- readRDS(bs1_files[i])
                  bs2 <- readRDS(bs2_files[i])
                  group_mn <- strsplit(basename(bs1_files[i]), "_bs1.rds")[[1]][1]
                  group_m <- strsplit(group_mn, "_vs_")[[1]][1]
                  group_n <- strsplit(group_mn, "_vs_")[[1]][2]

                  y_bs <- list()
                  y_bs[1] <- bs1
                  y_bs[2] <- bs2

                  this_m <- y_bs[[1]]
                  this_n <- y_bs[[2]]

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

                          log <- rbind(log,
                                paste0(group_mn, " bs_", j, ": ", nrow(adj_m)))
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
                      rownames(spectra_m) <- paste0(group_m, "_b", seqs)
                      rownames(spectra_n) <- paste0(group_n, "_b", seqs)
                      spectra_mn <- rbind(spectra_m, spectra_n)

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

                          log <- rbind(log,
                                paste0(group_mn, " bs_", j1, ": ", nrow(adj_m)))

                          m_j1 <- paste0(group_m, "_b", j1)

                          for (j2 in 1 : bs_len) {
                              adj_n <- unlist(this_n[[j2]])
                              adj_n[is.na(adj_n)] <- 0

                              n_j2 <- paste0(group_n, "_b", j2)

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
                      pm1 <- readRDS(pm1_files[i])
                      pm2 <- readRDS(pm2_files[i])

                      y_pm <- list()
                      y_pm[1] <- pm1
                      y_pm[2] <- pm2

                      this_mp <- y_pm[[1]]
                      this_np <- y_pm[[2]]

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

                              log <- rbind(log, paste0(group_mn, " pm_", k,
                                                ": ", nrow(adj_mp)))

                              if (nrow(adj_mp) < evk ||nrow(adj_np) < evk) {
                                  flag <- 1
                                  break
                              }

                              spectra_mp[k, ] <- get_spectra(adj_mp, k = evk)
                              spectra_np[k, ] <- get_spectra(adj_np, k = evk)
                          }

                          if (flag) next
                          seqs <- seq(1 : pm_len)
                          rownames(spectra_mp) <- paste0(group_m, "_p", seqs)
                          rownames(spectra_np) <- paste0(group_n, "_p", seqs)
                          spectra_mnp <- rbind(spectra_mp, spectra_np)
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

                              log <- rbind(log, paste0(group_mn, " pm_", k1,
                                                ": ", nrow(adj_mp)))

                              mp_k1 <- paste0(group_m, "_p", k1)

                              for (k2 in 1 : pm_len) {
                                  adj_np <- unlist(this_np[[k2]])
                                  adj_np[is.na(adj_np)] <- 0

                                  np_k2 <- paste0(group_n, "_p", k2)

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
                gc(reset = T)
              }
          write.table(log, paste0(out_dir, "log.txt"))

          y <- new("mina")
          y@dis_bs <- dis_bs
          y@dis_pm <- dis_pm

          ## get the stat table
          if (sig) {
              dis_stat <- get_stat(dis_bs, dis_pm)
          } else {
              dis_stat <- get_stat(dis_bs)
          }

          y@dis_stat <- dis_stat
          return(y)
        }
)
