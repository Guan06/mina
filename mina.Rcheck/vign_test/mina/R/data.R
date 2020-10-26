###############################################################################

#' Internal testing data of maize project, vegetative stage samples only,
#' including quantitative table (maize_asv.rds) and descriptive table
#' (maize_des.txt) for testing.
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
#' data(maize_asv)
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
#' data(maize_des)
#'
NULL

###############################################################################

#' Internal testing data of HMP project, including quantitative table (hmp_otu)
#' and descriptive table (hmp_des) for testing.
#'
#' @name data-hmp
#' @aliases hmp
#' @docType data
#' @keywords data
#' @examples
#' data(hmp)

###############################################################################
NA
###############################################################################

###############################################################################

#' OTU table of HMP project, data downloaded from
#' https://www.hmpdacc.org/hmp/HMQCP/
#'
#' @source HMP project.
#' @name hmp_otu
#' @format A matrix with samples in columns and OTUs in rows.
#' @examples
#' data(hmp_otu)
#'
NULL

#' Design file for HMP project, including 2711 samples in total.
#'
#' @source HMP project.
#' @name hmp_des
#' @format A data frame with columns:
#' \describe{
#'  \item{Sample_ID}{The unique ID of the microbial profiling sample.}
#'  \item{Sex}{The gender of the host human.}
#'  \item{Run_center}{The lab proccessing the sample sequencing.}
#'  \item{Subsite}{The subsite of body where samples were collected.}
#'  \item{Site}{The site of body where samples were collectec.}
#'  \item{Description}{The further details about the samples.}
#' }
#' @examples
#' data(hmp_des)
#'
NULL

###############################################################################

#' Subset of ASV table of maize project, ASVs appear in less than 100 samples
#' were filtered for later analysis.
#'
#' @source RECONSTRUCT project, maize microbiome part.
#' @name maize_asv2
#' @format A matrix with samples in columns and ASVs in rows. Unormalized table
#' including 1219 ASVs from 313 samples.
#' @examples
#' data(maize_asv2)
#'
NULL

#' Subset of design file of maize project, 313 samples are included.
#'
#' @source RECONSTRUCT project, maize microbiome part.
#' @name maize_des2
#' @format A data frame with columns:
#' \describe{
#'  \item{Sample_ID}{The unique ID of the microbial profiling sample.}
#'  \item{Host_genotype}{The genotype of the plant host maize.}
#'  \item{Compartment}{The compartment of the microbial sample comes from.}
#'  \item{Soil}{The soil of the sampled microbiome.}
#'  \item{Management}{The management of the soil where microbial sample from.}
#' }
#' @examples
#' data(maize_des2)
#'
NULL
