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
#' @return a (named) numeric vector. The length of the vector is equal to the
#' number of samples.
#'
#' @examples
#' example("PhIPData")
#' librarySize(phip_obj)
#'
#' ## Return an unnamed vector
#' librarySize(phip_obj, withDimnames = FALSE)
#' @export
librarySize <- function(object, ..., withDimnames = TRUE) {
    sums <- colSums(counts(object, ...))
    if (withDimnames) {
        sums
    } else {
        unname(sums)
    }
}

#' Proportion of sample reads
#'
#' This function calculates the proportion of total sample reads pulled by
#' each peptide.
#'
#' @param object \linkS4class{PhIPData} object
#' @param withDimnames logical; if true return a matrix with the same dimension
#' names as the original object.
#'
#' @return A (named) numeric matrix with the same dimensions as the function
#' input. Matrix values are between 0 and 1.
#'
#' @examples
#' example("PhIPData")
#' propReads(phip_obj)
#'
#' ## Return an unnamed matrix
#' propReads(phip_obj, withDimnames = FALSE)
#' @export
propReads <- function(object, withDimnames = TRUE) {
    n <- librarySize(object, withDimnames = FALSE)
    n_matrix <- matrix(rep(n, nrow(object)), nrow = nrow(object), byrow = TRUE)
    counts_matrix <- as.matrix(PhIPData::counts(object))
    prop_matrix <- counts_matrix / n_matrix
    if (withDimnames) {
        prop_matrix
    } else {
        unname(prop_matrix)
    }
}
