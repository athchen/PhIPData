#' @import methods SummarizedExperiment GenomicRanges IRanges
NULL

### ==============================================
### PhIPData class
### ==============================================
#' PhIPData - A class for PhIP-Seq experiment data
#'
#' @description The \code{PhIPData} class is used to manage the results from
#' phage-immunoprecipitation (PhIP-Seq) experiments. New \code{PhIPData} objects
#' can be created using the homonymous constructor.
#'
#' @return A \code{PhIPData} object
.PhIPData <- setClass("PhIPData", contains = "SummarizedExperiment")

### ==============================================
### PhIPData constructor
### ==============================================

#' Construct a \code{PhIPData} object.
#'
#' TODO Edit description. Empty objects are valid objects.
#' Missing peptideInfo and sampleInfo are also valid objects.
#'
#' @param counts a \code{matrix}, \code{data.frame}, or \code{\linkS4class{DataFrame}}
#' of integer read counts.
#' @param logfc a \code{matrix}, \code{data.frame}, or \code{\linkS4class{DataFrame}}
#' of log10 estimated fold changes.
#' @param prob a \code{matrix}, \code{data.frame}, or \code{\linkS4class{DataFrame}} of
#' probability values (p-values or posterior probabilities) for enrichment estimates.
#' @param peptideInfo a \code{data.frame} or \code{\linkS4class{DataFrame}} of peptide
#' information. Peptide identifiers (\code{pep_id}) can be specified via row names
#' of the \code{counts}, \code{logfc}, \code{prob}, or \code{peptideInfo} objects.
#' Valid construction of a \code{PhIPData} object only requires a \code{pep_id}, but
#' inclusion of other parameters, such as peptide positions or species, enables use
#' of other convenient functions for filtering peptides.
#' @param sampleInfo a \code{data.frame} or \code{\linkS4class{DataFrame}} of additional sample
#' information. Unique sample identifiers can be specific via column names of the
#' \code{counts}, \code{logfc}, \code{prob}, or \code{sampleInfo} objects.
#' @param .defaultNames vector of default names to use when the sample and peptide ID
#' disagree across the metadata and the \code{counts}, \code{logfc}, and \code{prob} matrices.
#' If \code{.defaultNames} is of length 1, the same source is used for both peptide and
#' sample identifiers. If \code{.defaultNames} is longer than 2, the first and second
#' elements correspond to the default names to use for peptides and samples, respectively.
#' Valid options are:
#' \itemize{
#'     \item{"info": }{names should be taken from the \code{SampleInfo} or \code{peptideInfo} objects.}
#'     \item{"counts": }{names should be taken from the row/column names of the \code{counts} object.}
#'     \item{"logfc": }{names should be taken from the row/column names of the \code{logfc} object.}
#'     \item{"prob": }{names should be taken from the row/column names of the \code{prob} object.}
#' }
#'
#' @export PhIPData
PhIPData <- function(counts = S4Vectors::DataFrame(),
                     logfc = S4Vectors::DataFrame(),
                     prob = S4Vectors::DataFrame(),
                     peptideInfo = S4Vectors::DataFrame(),
                     sampleInfo = S4Vectors::DataFrame(),
                     .defaultNames = "info", ...) {

  ## Variables defined for convenience
  arrays <- c("counts", "logfc", "prob")
  arrays_missing <- arrays[c(missing(counts), missing(logfc), missing(prob))]
  arrays_present <- arrays[!arrays %in% arrays_missing]
  num_missing <- length(arrays_missing)

  array_list <- list(counts = counts, logfc = logfc, prob = prob)
  .defaultNames <- if(length(.defaultNames) == 1) { rep(.defaultNames, 2) } else { .defaultNames }

  ## Check that input dimensions are matched.
  row_dims <- c(nrow(counts), nrow(logfc), nrow(prob), nrow(peptideInfo))
  col_dims <- c(ncol(counts), ncol(logfc), ncol(prob), nrow(sampleInfo))
  if(length(unique(row_dims[row_dims != 0])) > 1){
    stop("The number of samples is not consistent across inputs.")
  }
  if(length(unique(row_dims[row_dims != 0])) > 1){
    stop("The number of peptides is not consistent across inputs.")
  }

  ## Get peptide names
  peptide_names <- list(rownames(counts), rownames(logfc), rownames(prob), rownames(peptideInfo))
  peptide_names <- unique(peptide_names[sapply(peptide_names, length) != 0])
  peptide_warning <- "Peptide names are not identical across inputs. Using peptide names from "

  if (length(peptide_names) == 0) {
    peptide_names <- NULL
  } else if (length(peptide_names) == 1) {
    peptide_names <- peptide_names[[1]]
  } else if (.defaultNames[1] == "info") {
    warning(paste0(peptide_warning, "'peptideInfo'."))
    peptide_names <- rownames(peptideInfo)
  } else if(.defaultNames[1] == "counts") {
    warning(paste0(peptide_warning, "'counts'."))
    peptide_names <- colnames(counts)
  } else if (.defaultNames[1] == "logfc") {
    warning(paste0(peptide_warning, "'logfc'."))
    peptide_names <- colnames(logfc)
  } else if (.defaultNames[1] == "prob") {
    warning(paste0(peptide_warning, "'prob'."))
    peptide_names <- colnames(prob)
  } else {
    stop("Invalid '.defaultNames' supplied. Valid '.defaultNames' options are 'info', 'counts', 'logfc', or 'prob'.")
  }

  ## Get sample names
  sample_names <- list(colnames(counts), colnames(logfc), colnames(prob), rownames(sampleInfo))
  sample_names <- unique(sample_names[sapply(sample_names, length) != 0])
  sample_warning <- "Sample names are not identical across inputs. Using sample names from "

  if (length(sample_names) == 0) {
    sample_names <- NULL
  } else if (length(sample_names) == 1) {
    sample_names <- sample_names[[1]]
  } else if (.defaultNames[2] == "info") {
    warning(paste0(sample_warning, "'sampleInfo'."))
    sample_names <- rownames(sampleInfo)
  } else if(.defaultNames[2] == "counts") {
    warning(paste0(sample_warning, "'counts'."))
    sample_names <- rownames(counts)
  } else if (.defaultNames[2] == "logfc") {
    warning(paste0(sample_warning, "'logfc'"))
    sample_names <- rownames(logfc)
  } else if (.defaultNames[2] == "prob") {
    warning(paste0(sample_warning, "'prob'."))
    sample_names <- rownames(prob)
  } else {
    stop("Invalid '.defaultNames' supplied. Valid '.defaultNames' options are: 'info', 'counts', 'logfc', or 'prob'.")
  }

  ## Set missing arrays to DataFrames with dimensions and names corresponding to
  ## given matrices. All sample and peptide names are set to be identical
  ## even if they were mismatched in the inputs.
  if (num_missing == 3) {
    sapply(arrays_missing, function(array) array_list[[array]] <- S4Vectors::DataFrame())
  } else if (num_missing == 2) {
    empty_mat <- matrix(nrow = nrow(array_list[[arrays_present]]),
                        ncol = ncol(array_list[[arrays_present]]))
    rownames(empty_mat) <- peptide_names
    colnames(empty_mat) <- sample_names

    for(array in arrays_missing){
      array_list[[array]] <- S4Vectors::DataFrame(empty_mat)
    }

    array_list[[arrays_present]] <- S4Vectors::DataFrame(array_list[[arrays_present]])
    rownames(array_list[[arrays_present]]) <- peptide_names
    colnames(array_list[[arrays_present]]) <- sample_names

  } else if (num_missing == 1) {
    empty_mat <- matrix(nrow = nrow(array_list[[arrays_present[1]]]),
                        ncol = ncol(array_list[[arrays_present[1]]]))
    rownames(empty_mat) <- peptide_names
    colnames(empty_mat) <- sample_names

    array_list[[arrays_missing]] <- S4Vectors::DataFrame(empty_mat)
    for(array in arrays_present) {
      array_list[[array]] <- S4Vectors::DataFrame(array_list[[array]])
      rownames(array_list[[array]]) <- peptide_names
      colnames(array_list[[array]]) <- sample_names
    }
  }

  # Define peptide information
  if(missing(peptideInfo)){
    pep_meta <- S4Vectors::DataFrame(row.names = peptide_names)
    pep_start <- rep(0, length(peptide_names))
    pep_end <- rep(0, length(peptide_names))
  } else {
    pep_meta <- S4Vectors::DataFrame(peptideInfo[, !colnames(peptideInfo) %in% c("pos_start", "pos_end")])

    warnings <- c(no_start = (!"pos_start" %in% colnames(peptideInfo)),
                  no_end = (!"pos_end" %in% colnames(peptideInfo)))
    warnings["missing_start"] <- if(!warnings["no_start"]){
      any(is.na(peptideInfo[, "pos_start"])) } else { TRUE }
    warnings["missing_end"] <- if(!warnings["no_end"]){
      any(is.na(peptideInfo[, "pos_end"]))
    } else { TRUE }

    msg <- if (sum(warnings[c("no_start", "no_end")]) > 0) {
      paste0("Missing peptide ", paste(c("start", "end")[warnings[c("no_start", "no_end")]], collapse = " and "),
        " position information.")
    } else NULL

    msg <- if(any(warnings)) { paste0(msg, " Replacing missing values with 0.")} else { msg }

    warning(msg)

    if(!all(c("pos_start", "pos_end") %in% colnames(peptideInfo))){
      pep_start <- rep(0, length(peptide_names))
      pep_end <- rep(0, length(peptide_names))
    } else {
      pep_start <- if("pos_start" %in% colnames(peptideInfo)){
        replace(peptideInfo[, "pos_start"], is.na(peptideInfo[, "pos_start"]), 0)
        } else {
        rep(0, length(peptide_names))
        }

      pep_end <- if("pos_end" %in% colnames(peptideInfo)){
        replace(peptideInfo[, "pos_end"], is.na(peptideInfo[, "pos_end"]), 0)
      } else {
        rep(0, length(peptide_names))
      }

    }

  }

  row_info <- GenomicRanges::GRanges(seqnames = peptide_names,
                                     ranges = IRanges::IRanges(start = pep_start, end = pep_end))
  if(!S4Vectors::isEmpty(pep_meta)){ mcols(row_info) <- pep_meta }

  # Define sample information
  if(missing(sampleInfo)){
    sample_meta <- S4Vectors::DataFrame(row.names = sample_names)
  } else {
    sample_meta <- S4Vectors::DataFrame(sampleInfo)
  }

  # Make RangedSummarizedExperiment
  se_object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = array_list[["counts"]],
                  logfc = array_list[["logfc"]],
                  prob = array_list[["prob"]]))

  .PhIPData(se_object)
}


### ==============================================
### Getters
### ==============================================
setGeneric("counts", function(x) standardGeneric("counts"))
setGeneric("logfc", function(x) standardGeneric("logfc"))
setGeneric("prob", function(x) standardGeneric("prob"))
setGeneric("peptideInfo", function(x) standardGeneric("peptideInfo"))
setGeneric("sampleInfo", function(x) standardGeneric("sampleInfo"))

### ==============================================
### Setters
### ==============================================

### ==============================================
### Coercion methods
### ==============================================
