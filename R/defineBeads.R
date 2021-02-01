#' @export
getBeadsName <- function(){
  beads_name <- Sys.getenv("BEADS_NAME", "")

  if(beads_name == ""){ "beads" } else { beads_name }
}

#' @export
setBeadsName <- function(name){

  if(length(name) > 1){
    cli_alert_warning(paste0("Input has length larger than one. ",
                             "Using only the first element."))}
  if(typeof(name) != "character"){
    cli_alert_warning(paste0("Input is of type ", typeof(name), ". ",
                             "Coercing to character."))
  }

  name <- as.character(name[1])

  if(is.na(name)){
    stop("Beads cannot be specified via NA.")
  } else {
    Sys.setenv(BEADS_NAME = name)
  }
}
