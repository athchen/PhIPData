#' Peptide libraries
#'
#' @description PhIP-Seq experiments often use identical peptide libraries
#' different cohorts. These functions enable the user to conveniently reuse
#' tidied libraries.
#'
#' @details Each library is stored as a \linkS4class{DataFrame} in .rds file.
#' By default the libraries are stored in the \code{libraries} package folder.
#' \code{setLibraryPath} allows the user to change the location of this library,
#' if desired.
#'
#' While libraries can be directly saved in the specified format by the user,
#' it is highly recommended to use the \code{makeLibrary} function to
#' save a new library.
#'
#' @param path path to a folder containing a .rds files for each library
#' @param name name of the library
#' @param library a \code{matrix}, \code{data.frame}, or
#'     \linkS4class{DataFrame} with the peptide information for the
#'     specified library.
#'
#' @return \code{getLibraryPath} returns the path to the directory containing
#' .rds files for each library. \code{getLibrary} returns a
#' \linkS4class{DataFrame} corresponding to the peptide information for the
#' specified library.
#'
#' @examples
#' ## Get and set path to libraries folder
#' getLibraryPath()
#' \dontrun{
#' setLibraryPath("examplepath/")
#' }
#'
#' ## Create a new library
#' pep_meta <- data.frame(species = c(
#'     rep("human immunodeficiency virus", 3),
#'     rep("Epstein-Barr virus", 2)
#' ))
#' makeLibrary(pep_meta, "new_library")
#'
#' ## Use new library
#' counts_dat <- matrix(1:10, nrow = 5)
#' phip_obj <- PhIPData(
#'     counts = counts_dat,
#'     peptideInfo = getLibrary("new_library")
#' )
#'
#' ## Delete created library
#' file.remove(paste0(getLibraryPath(), "/new_library.rds"))
#' @name peptideLibraries
#' @include PhIPData-class.R
NULL

#' @describeIn peptideLibraries return the path to a folder containing the
#'     libraries.
#' @export
getLibraryPath <- function() {
    get("PHIP_LIBRARY_PATH", envir = pkg_env)
}

#' @describeIn peptideLibraries set the path to a folder containing the
#'     libraries.
#' @export
setLibraryPath <- function(path) {
    if (!is.character(path) | !dir.exists(path)) {
        stop("Invalid specified path.")
    } else {
        assign("PHIP_LIBRARY_PATH", normalizePath(path), envir = pkg_env)
        save(
            list = c("BEADS_NAME", "ALIAS_PATH", "PHIP_LIBRARY_PATH"),
            envir = pkg_env,
            file = system.file(package = "PhIPData", "extdata", "defaults.rda")
        )
    }
}

#' @describeIn peptideLibraries return a \linkS4class{DataFrame} with the
#' peptide information corresponding to the library.
#'
#' @export
getLibrary <- function(name) {
    path <- paste0(getLibraryPath(), "/", name, ".rds")
    readRDS(path)
}

#' @describeIn peptideLibraries create and store a \linkS4class{DataFrame} with
#' the specified peptide information.
#'
#' @export
#' @importFrom utils menu
makeLibrary <- function(library, name) {
    path <- paste0(getLibraryPath(), "/", name, ".rds")

    write <- if (file.exists(path)) {
        menu(c("Yes", "No"),
            title = paste0(
                "The library files already exists. ",
                "Do you want to overwrite the file?"
            )
        )
    } else {
        1
    }
    if (write == 1) {
        saveRDS(library, file = path)
    }
}
