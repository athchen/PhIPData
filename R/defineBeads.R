#' Defining how beads-only samples are encoded.
#'
#' @description \code{getBeadsName} and \code{setBeadsName} are two function to
#' get and set the string that encodes which samples are beads-only samples.
#' Information about beads-only samples are stored in the \code{groups} column
#' of \code{sampleInfo}.
#'
#' @details If \code{name} is of length greater than one, only the first element
#' of the vector is used. Non-character values of \code{name} are first coerced
#' into strings.
#'
#' @examples
#' ## Returns the default string, "beads"
#' getBeadsName()
#'
#' ## Not run since it changes defaults/user settings
#' \dontrun{
#' setBeadsName("beads-only")
#' }
#'
#' @name defineBeads
NULL

#' @describeIn defineBeads function that returns a string corresponding to how
#' beads-only samples are encoded.
#'
#' @return a string indicating how beads-only samples are encoded.
#' @export
getBeadsName <- function() {
    get("BEADS_NAME", envir = pkg_env)
}

#' @describeIn defineBeads function to set the string that indicates which
#' samples are beads-only samples in the \code{groups} column of
#' \code{sampleInfo}.
#'
#' @param name a string indicating how beads-only samples are encoded.
#' @export
#' @importFrom cli cli_alert_warning
setBeadsName <- function(name) {
    if (length(name) > 1) {
        cli::cli_alert_warning(paste0(
            "Input has length larger than one. ",
            "Using only the first element."
        ))
    }
    if (typeof(name) != "character") {
        cli::cli_alert_warning(paste0(
            "Input is of type ", typeof(name), ". ",
            "Coercing to character."
        ))
    }

    name <- as.character(name[1])

    if (is.na(name)) {
        stop("Beads cannot be specified via NA.")
    } else {
        assign("BEADS_NAME", name, envir = pkg_env)
        save(
            list = c("BEADS_NAME", "ALIAS_PATH", "PHIP_LIBRARY_PATH"),
            envir = pkg_env,
            file = system.file(package = "PhIPData", "extdata/defaults.rda")
        )
    }
}
