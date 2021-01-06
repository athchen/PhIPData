#' @include PhIPData-class.R

#' @export
getAliasPath <- function(){
  path <- Sys.getenv("ALIAS_PATH", "")

  if (path == ""){
    system.file(package = "PhIPData", "extdata/alias.rda")
  } else {
    path
  }
}

#' @export
setAliasPath <- function(path){
  if(!is.character(path) | !file.exists(path)){
    stop("Invalid specified file")
  } else if (!grepl("(rda|RData)", path)) {
    stop("Invalid file type.")
  } else {
    Sys.setenv(ALIAS_PATH = path)
  }
}


#' @export
getAlias <- function(virus){
  if(!virus %in% alias$alias){
    stop("Virus does not exist in alias database.")
  } else {
    alias$pattern[alias$alias == virus]
  }
}

#' @export
setAlias <- function(virus, pattern){

  if(virus %in% alias$alias){
    if(alias$pattern[alias$alias == virus] == pattern) {
      stop("Alias already exists in the alias database.")
    } else {
      alias$pattern[alias$alias == virus] <- pattern
    }
  } else {
    alias <- rbind(alias,
                   data.frame(alias = virus, pattern = pattern))
  }

  # IS THIS TOO FRAGILE?
  alias_loc <- getAliasPath()
  save(list = c("alias"), file = alias_loc)
  load(alias_loc, envir=parent.env(environment()))
}

#' @export
deleteAlias <- function(virus){

  if(!virus %in% alias$alias){
    stop("Virus does not exist in the alias database.")
  } else {
    virus_index <- which(alias$alias == virus)
    alias <- alias[-virus_index, ]
  }

  # IS THIS TOO FRAGILE?
  alias_loc <- getAliasPath()
  save(list = c("alias"), file = alias_loc)
  load(alias_loc, envir=parent.env(environment()))
}
