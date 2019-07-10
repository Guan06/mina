###############################################################################

#' Internal testing data of maize project, vegetative stage samples only,
#' including quantitative table (maize_asv.rds) and descriptive table
#' (maize_des.csv) for testing.
#'
#' @name data-maize
#' @aliases maize
#' @docType data
#' @keywords data
#' @examples
#' data(maize)

###############################################################################
NA
###############################################################################

###############################################################################

#' ASV table of maize project, vegetative stage samples only.
#'
#' @source RECONSTRUCT project, maize microbiome part.
#' @name maize_asv
#' @format A matrix with samples in columns and ASVs in rows. Unormalized table
#' including 12765 ASVs from 420 samples.
#' @examples
#' \dontrun{
#' maize_asv
#' }
#'
NULL

#' Design file of maize project, vegetative stage samples only, including 528
#' samples in total.
#'
#' @source RECONSTRUCT project, maize microbiome part.
#' @name maize_des
#' @format A data frame with columns:
#' \describe{
#'  \item{Sample_ID}{The unique ID of the microbial profiling sample.}
#'  \item{Host_genotype}{The genotype of the plant host maize.}
#'  \item{Compartment}{The compartment of the microbial sample comes from.}
#'  \item{Soil}{The soil of the sampled microbiome.}
#'  \item{Management}{The management of the soil where microbial sample from.}
#' }
#' @examples
#' \dontrun{
#' maize_des
#' }
#'
NULL
