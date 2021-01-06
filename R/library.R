#' @include PhIPData-class.R
#'
#'
#'

#' @export
getLibraryPath <- function(){
  path <- Sys.getenv("LIBRARY_PATH", "")

  if (path == ""){
    system.file(package = "PhIPData", "libraries")
  } else {
    path
  }

}

#' @export
setLibraryPath <- function(path){

  if(!is.character(path) | !dir.exists(path)){
    stop("Invalid specified path.")
  } else {
    Sys.setenv(LIBRARY_PATH = normalizePath(path))
  }

}

#' @export
getLibrary <- function(library){
  path <- paste0(getLibraryPath(), "/", library, ".rds")
  readRDS(path)
}

#' @export
makeLibrary <- function(library, name){
  path <- paste0(getLibraryPath(), "/", name, ".rds")
  saveRDS(library, file = path)
}
