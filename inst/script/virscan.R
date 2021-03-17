#' VirScan Library
#'
#' Below is the code to create a miniature peptide library to demontrate some
#' functionality of the PhIPData class. The full VirScan library, as described
#' in Xu et. al, consists of peptides that tile proteins from all human viruses.
#' The full library annotations can be found in the \code{pep_anno} object of the
#' \code{data_raw/hiv_data.rda} found
#' \href{https://github.com/athchen/ktsp_paper/}{here}.
#'
#' @seealso https://github.com/athchen/ktsp_paper/
#' @seealso Xu, G. J., Kula, T., Xu, Q., Li, M. Z., Vernon, S. D., Ndung'u, T.,
#'      Ruxrungtham, K., Sanchez, J., Brander, C., Chung, R. T.,
#'      O'Connor, K. C., Walker, B., Larman, H. B., & Elledge, S. J. (2015).
#'      Viral immunology. Comprehensive serological profiling of human
#'      populations using a synthetic human virome. Science
#'      (New York, N.Y.), 348(6239), aaa0698.
#'      https://doi.org/10.1126/science.aaa0698

## Download data from the link below.
## "https://github.com/athchen/ktsp_paper/raw/master/data_raw/hiv_data.rda"

library(readr)
library(dplyr)

load("hiv_data.rda")

set.seed(20210106)
pep_anno %>%
    group_by(species) %>%
    sample_n(ceiling(n()/100)) %>%
    write_tsv(file = "inst/extdata/virscan.tsv")
