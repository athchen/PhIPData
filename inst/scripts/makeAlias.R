#' Script to generate `alias.rda` in `inst/extdata`
#'
#' The `alias.rda` file contains a single dataframe with alias-pattern pairs.
#' Below is code to generate the initial alias database. For more information on
#' working with the alias database, see \code{?aliases}.
alias <- data.frame(alias = c("EBV", "HIV", "HPV"),
                    pattern = c("Epstein-Barr", "Human immunodeficiency virus",
                                "Human papillomavirus"))

save(list = "alias", file = "inst/extdata/alias.rda")
