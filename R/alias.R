#' Using aliases to subset virus data
#'
#' @description Rather than typing out full viruses names or repeating
#' regexpressions, users can use aliases as a convenient tool to subset
#' \code{PhIPData} objects by viral species.
#'
#' @details Aliases are saved to an rda file containing only a \code{data.frame}
#' with two columnes: \code{alias} and \code{pattern}. The \code{alias} column
#' contains the alias while the \code{pattern} column contains the corresponding
#' regexpression of interest.
#'
#' The location of the alias database is returned and defined by
#' \code{getAliasPath} and \code{setAliasPath}, respectively. By default
#' \code{getAliasPath} points to the \code{extdata} package folder.
#'
#' Once an alias is added to the database, it can always be accessed once the
#' package is loaded. It is recommended to use the functions \code{setAlias}
#' and \code{deleteAlias} to edit the alias database rather than modify the
#' .rda file itself. If an alias already exists in the database, \code{setAlias}
#' replaces the matched pattern.
#'
#' @param path path to \code{alias.rda}
#' @param virus character vector of the alias
#' @param pattern regexpression corresponding to the alias
#'
#' @examples
#' ## Get and set path to alias.rda
#' getAliasPath()
#' \dontrun{setAliasPath("examplepath/alias.rda")}
#'
#' ## Edit and modify aliases in the database
#' setAlias("test_virus", "test_pattern")
#' getAlias("test_virus")
#' setAlias("test_virus", "test_pattern2")
#' getAlias("test_virus")
#' deleteAlias("test_virus")
#'
#' ## Example of how to subset HIV using `getAlias`
#' ## Often, it is useful to set the `ignore.case` of `grep`/`grepl` to TRUE.
#' counts_dat <- matrix(1:10, nrow = 5)
#' peptide_meta <- data.frame(species = c(rep("human immunodeficiency virus", 3),
#'      rep("Epstein-Barr virus", 2)))
#'
#' phip_obj <- PhIPData(counts = counts_dat, peptideInfo = peptide_meta)
#' subset(phip_obj, grepl(getAlias("HIV"), species, ignore.case = TRUE))
#' @name aliases
#'
#' @include PhIPData-class.R
NULL

#' @describeIn aliases return the path to the .rda file of aliases.
#' @export
getAliasPath <- function(){
  path <- Sys.getenv("ALIAS_PATH", "")

  if (path == ""){
    system.file(package = "PhIPData", "extdata/alias.rda")
  } else {
    path
  }
}

#' @describeIn aliases set the path to the .rda file of aliases.
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

globalVariables("alias")
alias_env <- new.env(parent = emptyenv())
load(getAliasPath(), envir = alias_env)

#' @describeIn aliases return a regexpression corresponding to the alias.
#' @export
getAlias <- function(virus){
  if(!virus %in% alias_env$alias$alias){
    stop("Virus does not exist in the alias database.")
  } else {
    alias_env$alias$pattern[alias_env$alias$alias == virus]
  }
}

#' @describeIn aliases define/modify the regexpression for an alias.
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

#' @describeIn aliases remove an alias from the database.
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
