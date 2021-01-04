#' @include PhIPData-class.R
#'
#'
#'

set_library_path <- function(path){
  path <- Sys.getenv("LIBRARY_PATH", "")

  if (path == ""){
    stop("Library path does not exist. Please use {.fun setLibraryPath()}")
  } else {
    path
  }
}

get_library_path <- function(path){

}

use_library <- function(lib_name, path){

}
