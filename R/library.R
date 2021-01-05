#' @include PhIPData-class.R
#'
#'
#'

getLibraryPath <- function(){
  path <- Sys.getenv("LIBRARY_PATH", "")

  if (path == ""){
    stop("Library path does not exist.")
  } else {
    path
  }

}

setLibraryPath <- function(path){

  if(!is.character(path) | !dir.exists(path)){
    stop("Invalid specified path.")
  } else {
    Sys.setenv(LIBRARY_PATH = normalizePath(path))
  }

}

getLibrary <- function(library){
  path <- paste0(getLibraryPath(), "/", library, ".rds")
  readRDS(path)
}

makeLibrary <- function(library, name){
  path <- paste0(getLibraryPath(), "/", name, ".rds")
  saveRDS(library, file = path)
}
