###############################################################################

#' Visulization of network distance, average distances are used for tile plot.
#'
#' @import ggplot2
#' @import stringr
#' @param x An object of `mina` with slot `dis_stat` defined.
#' @param sig If `TRUE`, indicating significant distance with gold guild.
#' @param d The distance to be plotted, could be "BS" or "PM".
#' @return p The plotted figure.
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- bs_pm(maize, group = "Compartment")
#' maize <- net_dis(maize, method = "Jaccard")
#' p <- net_dis_plot(maize, d = "BS")
#' @rdname net_dis_plot
#' @export

setMethod("net_dis_plot", signature("mina", "ANY"),
          function(x, d = "BS", sig = TRUE){
              
              dis <- dis_stat(x)
              dis$Sig[dis$p < 0.05] <- "sig"
              dis$Sig[dis$p >=0.05] <- "non-sig"

              if(d == "BS") {
                  dis <- dis[, c("Compare", "Distance_Mean", "Sig")]
              }else{
                  dis <- dis[, c("Compare", "Distance_PM_Mean", "Sig")]
              }
              colnames(dis)[2] <- "Distance"
    
              if(sig) {
                  c_sig <- c("sig" = "gold",
                             "non-sig" = "gray")
              }

              if (max(dis$Distance) <= 1) {di <- 2} else{di <- 0}

              dis$Group1 <- str_split_fixed(dis$Compare, "_", 2)[, 1]
              dis$Group2 <- str_split_fixed(dis$Compare, "_", 2)[, 2]

              p <- ggplot(dis, aes(Group1, Group2)) +
                geom_tile(aes(fill = Distance, color = Sig), size = 1.2) +
                geom_text(aes(label = round(Distance, di), color = Sig), size = 4) +
                scale_color_manual(values = c_sig) +
                theme(panel.background = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    axis.ticks = element_line(color = "black"),
                    axis.text.x = element_text(colour = "black", angle = 90,
                                     size = 8, hjust = 1),
                    axis.text.y = element_text(size = 8),
                    aspect.ratio=1,
                    text = element_text(family="sans"))
          }
)
