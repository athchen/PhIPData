#' Peptide libraries
#'
#' @description PhIP-Seq experiments often use identical peptide libraries
#' different cohorts. These functions enable the user to conveniently reuse
#' tidied libraries.
#'
#' @details Each library is stored as a \linkS4class{DataFrame} in .rds file.
#' New libraries can be stored for future use with the \code{makeLibrary}
#' function.
#'
#' @param name name of the library
#' @param library a \code{matrix}, \code{data.frame}, or
#'     \linkS4class{DataFrame} with the peptide information for the
#'     specified library.
#'
#' @return \code{getLibrary} returns a \linkS4class{DataFrame} corresponding to
#' the peptide information for the specified library.
#'
#' @examples
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
#' ## List libraries
#' listLibrary()
#'
#' ## Delete created library
#' removeLibrary("new_library")
#' @name peptideLibraries
#' @include PhIPData-class.R
NULL

#' @describeIn peptideLibraries return a \linkS4class{DataFrame} with the
#' peptide information corresponding to the library.
#'
#' @export
getLibrary <- function(name) {
    path <- BiocFileCache::bfcquery(
        pkg_env$beer_cache,
        paste0("library_", name)
    )$rpath
    readRDS(path)
}

#' @describeIn peptideLibraries create and store a \linkS4class{DataFrame} with
#' the specified peptide information.
#'
#' @export
makeLibrary <- function(library, name) {
    ## Check if library exists
    query_out <- BiocFileCache::bfcquery(
        pkg_env$beer_cache,
        paste0("library_", name)
    )$rpath

    ## If the library exists, confirm that the file should be overwritten
    write <- if (length(query_out) != 0) {
        utils::menu(c("Yes", "No"),
            title = paste0(
                "The library files already exists. ",
                "Do you want to overwrite the file?"
            )
        )
    } else {
        1
    }

    if (write == 1) {
        # Define library path
        library_path <- ifelse(
            length(query_out) == 0,
            BiocFileCache::bfcnew(
                pkg_env$beer_cache,
                paste0("library_", name),
                ext = ".rds"
            ),
            query_out
        )
        saveRDS(library, file = library_path)
    }
}

#' @describeIn peptideLibraries delete stored libraries
#'
#' @export
removeLibrary <- function(name) {

    ## Check if library exists
    query_out <- BiocFileCache::bfcquery(
        pkg_env$beer_cache,
        paste0("library_", name)
    )$rid

    if (length(query_out) == 0) {
        cli::cli_alert_warning("Library not found. No libraries removed.")
    } else {
        BiocFileCache::bfcremove(pkg_env$beer_cache, query_out)
    }
}


#' @describeIn peptideLibraries list all available libraries
#'
#' @export
listLibrary <- function() {

    ## Check if library exists
    BiocFileCache::bfcquery(
        pkg_env$beer_cache,
        "library_"
    )$rname
}
