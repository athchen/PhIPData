#' @include PhIPData-class.R
NULL

#' Subset beads-only samples
#'
#' Function to subset PhIP-seq data for beads-only samples.
#'
#' @param object \linkS4class{PhIPData} object
#' @return a \linkS4class{PhIPData} object.
#'
#' @examples
#' example("PhIPData")
#' subsetBeads(phip_obj)
#' @export
subsetBeads <- function(object) {
    object[, object$group == getBeadsName()]
}
