###############################################################################

#' Visulization of components distance / dissimilarity in k dimension.
#'
#' @import ggplot2
#' @param x An object of `mina` with list `dmr` defined.
#' @param match The column name of the components IDs in `des` with exactly the
#' same as rownames in x.
#' @param d1 The dimension be visualized in x-axis, default `1`.
#' @param d2 The dimension be visualized in y-axis, default `2`.
#' @param color The column name in `des` to be used for different color groups.
#' @param shape The column name in `des` to be used for different shape groups,
#' default `NULL`.
#' @param ... Additional parameters.
#' @return p The plotted figure.
#' @examples
#' maize <- new("mina", tab = maize_asv, des = maize_des)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' p1a <- com_plot(maize, match = "Sample_ID", color = "Compartment")
#' p1b <- com_plot(maize, match = "Sample_ID", d1 = 3, d2 = 4,
#' color = "Compartment")
#' p2a <- com_plot(maize, match = "Sample_ID", color = "Host_genotype")
#' p2b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 3, color =
#' "Host_genotype")
#' p3a <- com_plot(maize, match = "Sample_ID", color = "Compartment", shape =
#' "Soil")
#' p3b <- com_plot(maize, match = "Sample_ID", d1 = 1, d2 = 4, color =
#' "Compartment", shape = "Soil")
#' @return The PCoA plot.
#' @rdname com_plot-mina
#' @exportMethod com_plot

setMethod("com_plot", signature("mina", "ANY"),
          function(x, match, d1 = 1, d2 = 2, color, shape = NULL, ...) {
              stop("Mush specify a `match` column.")
          }
)
###############################################################################

#' @rdname com_plot-mina
#' @exportMethod com_plot

setMethod("com_plot",
          signature("mina", "character"),
          function(x, match, d1 = 1, d2 = 2, color, shape = NULL, ...) {
              stop("Mush specify a column for `color`.")
          }
)

###############################################################################

#' @rdname com_plot-mina
#' @exportMethod com_plot

setMethod("com_plot",
          signature("mina", "character"),
          function(x, match, d1 = 1, d2 = 2, color, shape = NULL, ...) {
              stopifnot(
                    is.character(c(match, color)),
                    is.numeric(c(d1, d2))
              )
              p <- pcoa_plot(.dmr(x), des(x), match = match, d1 = d1, d2 = d2,
                             color = color, shape = shape)
              return(p)
          }
)

###############################################################################

#' Visulization of components distance / dissimilarity in k dimension.
#'
#' @import ggplot2
#' @param x A list generated by `dmr`.
#' @param des The corresponding descriptive table.
#' @param match The column name of the components IDs in `des` with exactly the
#' same as rownames in x.
#' @param d1 The dimension be visualized in x-axis, default `1`.
#' @param d2 The dimension be visualized in y-axis, default `2`.
#' @param color The column name in `des` to be used for different color groups.
#' @param shape The column name in `des` to be used for different shape groups,
#' default `NULL`.
#' @param ... Additional parameters.
#' @return p The plotted PCoA.
#' @rdname pcoa_plot
#' @examples
#' maize <- new("mina", tab = maize_asv2, des = maize_des2)
#' maize <- norm_tab(maize, method = "raref")
#' maize <- fit_tabs(maize)
#' maize <- com_dis(maize, method = "bray")
#' maize <- dmr(maize)
#' asv_dmr <- .dmr(maize)
#' des <- des(maize)
#' p1a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment")
#' p1b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 3, d2 = 4, color =
#' "Compartment")
#' p2a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Host_genotype")
#' p2b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 3, color =
#' "Host_genotype")
#' p3a <- pcoa_plot(asv_dmr, des, match = "Sample_ID", color = "Compartment",
#' shape = "Soil")
#' p3b <- pcoa_plot(asv_dmr, des, match = "Sample_ID", d1 = 1, d2 = 4, color =
#' "Compartment", shape = "Soil")
#' @exportMethod pcoa_plot

setMethod("pcoa_plot", signature("list", "data.frame", "character"),
            function(x, des, match, d1 = 1, d2 = 2,
                     color, shape = NULL, ...) {
    
                stopifnot(
                          is.character(c(match, color)),
                          is.numeric(c(d1, d2))
                )

    points <- x$points
    if (nrow(points) != nrow(des)) {
        stop("Component number in `dmr` and `des` are different!")
    }
    if (all(sort(rownames(points)) != sort(des[[match]]))) {
        stop("Component IDs in `dmr` and `des` are different!")
    }

    # finish checking and start plotting
    points <- as.data.frame(points)
    colnames(points)[c(d1, d2)] <- c("x", "y")

    des <- des[match(rownames(points), des[[match]]), ]
    points <- cbind(points, des[[color]])
    colnames(points)[ncol(points)] <- color
    if (! is.null(shape)) {
        points <- cbind(points, des[[shape]])
        colnames(points)[ncol(points)] <- shape
    }

    eig <- x$eig
    eig[eig < 0] <- 0
    eig1 <- 100 * eig[d1] / sum(eig)
    eig2 <- 100 * eig[d2] / sum(eig)

    p <- ggplot(points, aes_string(x = "x", y = "y", color = color,
                                   shape = shape)) +
                geom_point(size = 1.2, alpha = 0.8) +
                coord_fixed(ratio = 1) +
                labs (x = paste0("PCo", d1, " (",
                                 format(eig1, digits = 4), "%)"),
                      y = paste0("PCo", d2, " (",
                                 format(eig2, digits = 4), "%)")) +
                theme(panel.background = element_blank(),
                      panel.grid = element_blank(),
                      axis.line.x = element_line(color = "black"),
                      axis.line.y = element_line(color = "black"),
                      axis.ticks = element_line(color = "black"),
                      axis.text = element_text(colour = "black", size = 12),
                      legend.position = "right",
                      legend.background = element_blank(),
                      legend.key = element_blank(),
                      text = element_text(family="sans")
                )
    return(p)
})
