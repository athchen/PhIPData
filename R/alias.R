#' Using aliases to subset virus data
#'
#' @description Rather than typing out full viruses names or repeating
#' regexpressions, users can use aliases as a convenient tool to subset
#' \code{PhIPData} objects by viral species.
#'
#' @details Aliases are saved to an rda file containing only a \code{data.frame}
#' with two columns: \code{alias} and \code{pattern}. The \code{alias} column
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
#' replaces the matched pattern. If an alias does not exist in the database,
#' \code{getAlias} returns \code{NA_character_}.
#'
#' @param path path to \code{alias.rda}
#' @param virus character vector of the alias
#' @param pattern character vector of regexpressions corresponding to the alias
#'
#' @return \code{getAliasPath()} returns the path to the alias database.
#' \code{getAlias()} returns a vector of regexpressions corresponding to
#' queried inputs. The returned vector is the same length as the input vector.
#' Queries that do not exist in the database return \code{NA_character_}.
#'
#' @examples
#' ## Get and set path to alias.rda
#' getAliasPath()
#' \dontrun{
#' setAliasPath("examplepath/alias.rda")
#' }
#'
#' ## Edit and modify aliases in the database
#' setAlias("test_virus", "test_pattern")
#' getAlias("test_virus")
#' setAlias("test_virus", "test_pattern2")
#' getAlias("test_virus")
#' deleteAlias("test_virus")
#'
#' ## Edit and modify multiple aliases at once.
#' setAlias(c("virus_1", "virus_2"), c("pattern_1", "pattern_2"))
#' getAlias(c("virus_1", "virus_2"))
#' deleteAlias(c("virus_1", "virus_2"))
#'
#' ## Example of how to subset HIV using `getAlias`
#' ## Often, it is useful to set the `ignore.case` of `grep`/`grepl` to TRUE.
#' counts_dat <- matrix(1:10, nrow = 5)
#' peptide_meta <- data.frame(species = c(rep("Epstein-Barr virus", 2),
#'     rep("human immunodeficiency virus", 3)))
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

.getOneAlias <- function(virus){
  if(!virus %in% alias_env$alias$alias){ NA_character_ } else {
    alias_env$alias$pattern[alias_env$alias$alias == virus]
  }
}

#' @describeIn aliases return a regexpression corresponding to the alias.
#' @export
getAlias <- function(virus){
  vapply(virus, .getOneAlias, character(1), USE.NAMES = FALSE)
}


#' @describeIn aliases define/modify the regexpression for an alias.
#' @export
setAlias <- function(virus, pattern){
  if(length(virus)!= length(pattern)){
    stop("Input vector lengths are unequal.")
  }

  current_alias <- alias_env$alias

  ## Look at whether any viruses need to be added or changed
  ##    new_viruses: viruses to be added
  ##    exist_viruses: viruses that exist in the database
  ##        (may have the same patterns)
  ##    replace_viruses: subset of exist viruses that need to have the pattern
  ##        changed
  new_viruses <- virus[!virus %in% current_alias$alias]
  exist_viruses <- setdiff(virus, new_viruses)
  current_pattern <- vapply(exist_viruses, function(x) {
    pattern <- current_alias$pattern[current_alias$alias == x]
    if(length(pattern) == 0) NA else pattern
    }, character(1))
  replace_viruses <- exist_viruses[current_pattern !=
                                     pattern[virus == exist_viruses]]
  n_replace <- length(replace_viruses)
  if(n_replace > 0){
    cli::cli_alert_info("Replacing pattern{?s} for {n_replace} alias{?es}.")
  }

  if(sum(length(new_viruses), n_replace) ==0){
    cli::cli_alert_info("No new alias-pattern combinations added.")
  } else {
    ## Add new aliases
    new_alias <- rbind(current_alias,
                       data.frame(alias = new_viruses,
                                  pattern = pattern[virus == new_viruses]))
    ## Replace patterns
    for(replacement in replace_viruses){
      new_alias$pattern[new_alias$alias == replacement] <-
        pattern[virus == replacement]
    }

    ## Change environment
    alias_env$alias <- new_alias

    ## Save to where the environment is loaded
    alias_env$alias_loc <- getAliasPath()
    save(alias,
         envir = alias_env,
         file = alias_env$alias_loc)
  }
}

.deleteOneAlias <- function(virus){
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

#' @describeIn aliases remove an alias from the database.
#' @export
deleteAlias <- function(virus){
  for(i in virus) .deleteOneAlias(i)
}
