#' @include PhIPData-class.R alias.R

#' Subset a PhIPData object by virus species.
#'
#' @param object a \code{PhIPData} object
#' @param virus regex indicating the virus of interest
#' @return a \code{PhIPData} object
getPeptides <- function(object, virus) {
  if(!"species" %in% colnames(mcols(peptideInfo(object)))){
    stop("Peptide metadata does not contain `species` information.")
  } else {
    peptide_ind <- which(grepl(virus, mcols(peptideInfo(object))[, "species"]))

    PhIPData(counts = counts(object)[peptide_ind, ],
             logfc = logfc(object)[peptide_ind, ],
             prob = logfc(object)[peptide_ind, ],
             peptideInfo = peptideInfo(object)[peptide_ind, ],
             sampleInfo = sampleInfo(object),
             metadata = metadata(object))
  }
}
