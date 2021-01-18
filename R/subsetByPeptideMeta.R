#' @include PhIPData-class.R

#' @title Subset a PhIPData object by virus species.
#'
#' @param object a \code{PhIPData} object
#' @param virus regex indicating the virus of interest
#' @return a \code{PhIPData} object
#' @export
#' @importFrom dplyr as_tibble mutate filter pull
subsetByPeptideMeta <- function(x, ...) {

  pep_meta <- peptideInfo(x) %>%
    mcols() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(row_number = 1:dplyr::n())

  pep_ind <- dplyr::filter(pep_meta, ...) %>%
    dplyr::pull(row_number)

  x[pep_ind, ]
}
