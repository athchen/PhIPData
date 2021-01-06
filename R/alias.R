#' @include PhIPData-class.R

getAliasPath <- function(){
  path <- Sys.getenv("ALIAS_PATH", "")

  if (path == ""){
    system.file(package = "PhIPData", "R/sysdata.rda")
  } else {
    path
  }
}

setAliasPath <- function(path){
  if(!is.character(path)){
    stop("Invalid specified path.")
  } else {
    Sys.setenv(ALIAS_PATH = path)
  }
}

getAlias <- function(virus){
  if(!virus %in% alias$alias){
    stop("Virus does not exist in alias database.")
  } else {
    alias$pattern[alias$alias == virus]
  }
}

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
