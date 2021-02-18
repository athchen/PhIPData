#' @include PhIPData-class.R
NULL

#' Calculate total read counts for each sample.
#'
#' This function is a wrapper function for \link{colSums} on the \code{counts}
#' assay.
#'
#' @param object \linkS4class{PhIPData} object
#' @param ... arguments passed to \link{colSums}
#' @param withDimnames logical; if true, the vector names are the sample names;
#' otherwise the vector is unnamed.
#'
#' @return a named numeric vector. The length of the vector is equal to the
#' number of samples.
#'
#' @examples
#' example("PhIPData")
#' librarySize(phip_obj)
#'
#' ## Return an unnamed vector
#' librarySize(phip_obj, withDimnames = FALSE)
#'
#' @export
librarySize <- function(object, ..., withDimnames = TRUE){
  sums <- colSums(counts(object,...))
  if(withDimnames){ sums } else { unname(sums) }
}
