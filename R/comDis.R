################################################################################
#' List of distance method keys supported in \code{\link{comDis}}
#'
#' Add methods can only be matched exactly.

#' @format A list of character vectors.
#'
#' \describe{
#'   \item{\code{tina}}{\code{\link[phyloseq]{UniFrac}}}
#'   \item{\code{manhattan}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{euclidean}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{canberra}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{bray}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{kulczynski}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{jaccard}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{gower}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{altGower}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{morisita}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{horn}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{mountford}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{raup}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{binomial}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{chao}}{\code{\link[vegan]{vegdist}}}
#'   \item{\code{cao}}{\code{\link[vegan]{vegdist}}}
#' }
#'
#' @export
#'
#' @examples
#' distanceMethodList
distanceMethodList <- list(
  # The methods supported by vegan::vegdist function.
  vegdist    = c("manhattan", "euclidean", "canberra", "bray",
                 "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",
                 "mountford", "raup" , "binomial", "chao", "cao"),

  TINA       = "tina",
  designdist = "ANY"
)
