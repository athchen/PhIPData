#' @include PhIPData-class.R
#'
#'
#'

#' @export
getLibraryPath <- function(){
  path <- Sys.getenv("PHIP_LIBRARY_PATH", "")

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
    Sys.setenv(PHIP_LIBRARY_PATH = normalizePath(path))
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

  write <- if(file.exists(path)){
    menu(c("Yes", "No"),
         title = "The library files already exists. Do you want to overwrite the file?")
  } else { 1 }
  if(write == 1) { saveRDS(library, file = path) }
}
