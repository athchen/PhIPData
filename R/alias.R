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
    stop("Invalid alias file location.")
  } else if (!grepl("(rda|RData)", path)) {
    stop("Invalid file type.")
  } else {
    Sys.setenv(ALIAS_PATH = path)
  }
}


alias_env <- new.env(parent = emptyenv())
load(getAliasPath(), envir = alias_env)

#' @export
getAlias <- function(virus){
  if(!virus %in% alias_env$alias$alias){
    stop("Virus does not exist in the alias database.")
  } else {
    alias_env$alias$pattern[alias_env$alias$alias == virus]
  }
}

#' @export
setAlias <- function(virus, pattern){

  if(virus %in% alias_env$alias$alias){

    if(alias_env$alias$pattern[alias_env$alias$alias == virus] == pattern) {
      stop("Alias already exists in the alias database.")
    } else {
      alias_env$alias$pattern[alias_env$alias$alias == virus] <- pattern
    }

  } else {
    alias_env$alias <- rbind(alias_env$alias,
                             data.frame(alias = virus, pattern = pattern))
  }

  alias_env$alias_loc <- getAliasPath()
  save(alias,
       envir = alias_env,
       file = alias_env$alias_loc)
}

#' @export
deleteAlias <- function(virus){

  if(!virus %in% alias_env$alias$alias){
    stop("Virus does not exist in the alias database.")
  } else {
    virus_index <- which(alias_env$alias$alias == virus)
    alias_env$alias <- alias_env$alias[-virus_index, ]
  }

  alias_env$alias_loc <- getAliasPath()
  save(alias,
       envir = alias_env,
       file = alias_env$alias_loc)
}
