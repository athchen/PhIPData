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
#' Missing peptideInfo and sampleInfo are also valid objects, empty assays can also be initialized
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
  assays <- c("counts", "logfc", "prob")
  assays_missing <- assays[c(missing(counts), missing(logfc), missing(prob))]
  assays_present <- assays[!assays %in% assays_missing]
  num_missing <- length(assays_missing)

  assay_list <- list(counts = counts, logfc = logfc, prob = prob)
  .defaultNames <- if(length(.defaultNames) == 1) { rep(.defaultNames, 2) } else { .defaultNames }

  ## Check that input dimensions are matched.
  dims <- .checkDims(counts, logfc, prob, peptideInfo, sampleInfo)
  if(is.character(dims)) { stop(dims) }

  # ## Check that both peptide and sample information are present
  # if(num_missing == 3 & sum(c(missing(sampleInfo), missing(peptideInfo))) == 1) {
  #   stop("Cannot create empty PhIPData object with only one of sampleInfo or peptideInfo.")
  # }

  ## Get peptide names
  peptide_names <- .getPeptideNames(counts, logfc, prob, peptideInfo, .defaultNames[1])

  ## Get sample names
  sample_names <- .getSampleNames(counts, logfc, prob, sampleInfo, .defaultNames[2])

  ## Set missing assays to DataFrames with dimensions and names corresponding to
  ## given matrices. All sample and peptide names are set to be identical
  ## even if they were mismatched in the inputs.
  assay_list <- .tidyAssays(assay_list, assays_missing, peptide_names, sample_names)

  ## Define peptide information
  tidied_pepInfo <- .tidyPeptideInfo(peptideInfo, peptide_names)
  pep_meta <- tidied_pepInfo[["pep_meta"]]
  pep_start <- tidied_pepInfo[["pep_start"]]
  pep_end <- tidied_pepInfo[["pep_end"]]

  row_info <- GenomicRanges::GRanges(seqnames = peptide_names,
                                     ranges = IRanges::IRanges(start = pep_start, end = pep_end))
  if(!S4Vectors::isEmpty(pep_meta)){ mcols(row_info) <- pep_meta }

  ## Define sample information
  if(missing(sampleInfo)){
    sample_meta <- S4Vectors::DataFrame(row.names = sample_names)
  } else {
    sample_meta <- S4Vectors::DataFrame(sampleInfo, row.names = sample_names)
  }

  ## Make RangedSummarizedExperiment
  se_object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = assay_list[["counts"]],
                  logfc = assay_list[["logfc"]],
                  prob = assay_list[["prob"]]))

  .PhIPData(se_object)
}


.checkDims <- function(counts, logfc, prob, peptideInfo, sampleInfo){

  row_dims <- c(nrow(counts), nrow(logfc), nrow(prob), nrow(peptideInfo))
  col_dims <- c(ncol(counts), ncol(logfc), ncol(prob), nrow(sampleInfo))

  fixed_row <- if(length(unique(row_dims)) == 1) { unique(row_dims) } else { unique(row_dims[row_dims != 0]) }
  fixed_col <- if(length(unique(col_dims)) == 1) { unique(col_dims) } else { unique(col_dims[col_dims != 0]) }

  match <- c(length(fixed_row) == 1, length(fixed_col) == 1)

  error <- paste0("The number of ",
                  paste0(c("samples", "peptides")[! match], collapse = " and "),
                  " differs across inputs.")

  if(sum(match) < 2){ error } else { c(fixed_row, fixed_col) }
}

.getPeptideNames <- function(counts, logfc, prob, peptideInfo, default){

  peptide_names <- list(rownames(counts), rownames(logfc), rownames(prob), rownames(peptideInfo))
  peptide_names <- unique(peptide_names[sapply(peptide_names, length) != 0])

  peptide_warning <- "Peptide names are not identical across inputs. Using peptide names from "
  default_error <- "Invalid '.defaultNames' supplied. Valid '.defaultNames' options are 'info', 'counts', 'logfc', or 'prob'."

  if (length(peptide_names) == 0) {
    peptide_names <- NULL
  } else if (length(peptide_names) == 1) {
    peptide_names <- peptide_names[[1]]
  } else if (default == "info") {
    warning(paste0(peptide_warning, "'peptideInfo'."))
    peptide_names <- rownames(peptideInfo)
  } else if(default == "counts") {
    warning(paste0(peptide_warning, "'counts'."))
    peptide_names <- colnames(counts)
  } else if (default== "logfc") {
    warning(paste0(peptide_warning, "'logfc'."))
    peptide_names <- colnames(logfc)
  } else if (default == "prob") {
    warning(paste0(peptide_warning, "'prob'."))
    peptide_names <- colnames(prob)
  } else { stop(default_error) }

  peptide_names
}

.getSampleNames <- function(counts, logfc, prob, sampleInfo, default){

  sample_names <- list(colnames(counts), colnames(logfc), colnames(prob), rownames(sampleInfo))
  sample_names <- unique(sample_names[sapply(sample_names, length) != 0])

  sample_warning <- "Sample names are not identical across inputs. Using sample names from "
  default_error <- "Invalid '.defaultNames' supplied. Valid '.defaultNames' options are: 'info', 'counts', 'logfc', or 'prob'."

  if (length(sample_names) == 0) {
    sample_names <- NULL
  } else if (length(sample_names) == 1) {
    sample_names <- sample_names[[1]]
  } else if (default == "info") {
    warning(paste0(sample_warning, "'sampleInfo'."))
    sample_names <- rownames(sampleInfo)
  } else if(default == "counts") {
    warning(paste0(sample_warning, "'counts'."))
    sample_names <- rownames(counts)
  } else if (default == "logfc") {
    warning(paste0(sample_warning, "'logfc'"))
    sample_names <- rownames(logfc)
  } else if (default == "prob") {
    warning(paste0(sample_warning, "'prob'."))
    sample_names <- rownames(prob)
  } else { stop(default_error) }

  sample_names
}

.tidyAssays <- function(assay_list, assays_missing, peptide_names, sample_names){

  assays <- c("counts", "logfc", "prob")
  assays_present <- assays[!assays %in% assays_missing]
  num_missing <- length(assays_missing)

  if (num_missing == 3) {
    sapply(assays_missing, function(assay) assay_list[[assay]] <- S4Vectors::DataFrame())
  } else if (num_missing == 2) {
    empty_mat <- matrix(nrow = nrow(assay_list[[assays_present]]),
                        ncol = ncol(assay_list[[assays_present]]))
    rownames(empty_mat) <- peptide_names
    colnames(empty_mat) <- sample_names

    for(assay in assays_missing){
      assay_list[[assay]] <- S4Vectors::DataFrame(empty_mat)
    }

    assay_list[[assays_present]] <- S4Vectors::DataFrame(assay_list[[assays_present]])
    rownames(assay_list[[assays_present]]) <- peptide_names
    colnames(assay_list[[assays_present]]) <- sample_names

  } else if (num_missing == 1) {
    empty_mat <- matrix(nrow = nrow(assay_list[[assays_present[1]]]),
                        ncol = ncol(assay_list[[assays_present[1]]]))
    rownames(empty_mat) <- peptide_names
    colnames(empty_mat) <- sample_names

    assay_list[[assays_missing]] <- S4Vectors::DataFrame(empty_mat)
    for(assay in assays_present) {
      assay_list[[assay]] <- S4Vectors::DataFrame(assay_list[[assay]])
      rownames(assay_list[[assay]]) <- peptide_names
      colnames(assay_list[[assay]]) <- sample_names
    }
  }

  assay_list
}

.tidyPeptideInfo <- function(peptideInfo, peptide_names){
  if(all(dim(peptideInfo) == 0)){
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

  list(pep_meta = pep_meta,
       pep_start = pep_start,
       pep_end = pep_end)
}

### ==============================================
### Validity
### ==============================================

## 1. counts cannot have negative entries
.checkAssays <- function(object){
  if(any(counts(object) < 0)) { "'counts' cannot have negative entries." } else NULL
}

## 2. Sample and peptide dimensions must be the same within the object.
.checkObjectDim <- function(object){
  dim_check <- .checkDims(counts(object), logfc(object), prob(object), peptideInfo(object), sampleInfo(object))
  if(is.character(dim_check)) { dim_check } else NULL
}

## 3. sample and peptide names must be identical across all assays and annotation information.
.checkNames <- function(object){

  sample_names <- c(colnames(counts(object)), colnames(logfc(object)),
                    colnames(prob(object)), rownames(sampleInfo(object)))
  peptide_names <- c(rownames(counts(object)), rownames(logfc(object)),
                     rownames(prob(object)), rownames(peptideInfo(object)))

  match <- c(length(unique(sample_names)) == 1, length(unique(peptide_names)) == 1)

  error <- paste0("Names do not match across ",
                  paste0(c("samples", "peptides")[!match], collapse = " and "), ".")

  if(sum(match < 2)) { error } else { NULL }

}

.validPhIPData <- function(object){
  if(!isEmpty(object)){
    c(.checkAssays(object), .checkObjectDim(object), .checkNames(object))
  }
}

S4Vectors::setValidity2("PhIPData", .validPhIPData)

### ==============================================
### Getters
### ==============================================

setGeneric("counts", function(x) standardGeneric("counts"))
setMethod("counts", "PhIPData", function(x) assays(x)[["counts"]])

setGeneric("logfc", function(x) standardGeneric("logfc"))
setMethod("logfc", "PhIPData", function(x) assays(x)[["logfc"]])

setGeneric("prob", function(x) standardGeneric("prob"))
setMethod("prob", "PhIPData", function(x) assays(x)[["prob"]])

setGeneric("peptideInfo", function(x) standardGeneric("peptideInfo"))
setMethod("peptideInfo", "PhIPData", function(x) rowRanges(x))

setGeneric("sampleInfo", function(x) standardGeneric("sampleInfo"))
setMethod("sampleInfo", "PhIPData", function(x) colData(x))

### ==============================================
### Setters
### ==============================================

setGeneric("counts<-", function(object, value) standardGeneric("counts<-"))
setReplaceMethod("counts", "PhIPData", function(object, value) {
  .replaceAssay(object, "counts", value)
})

setGeneric("logfc<-", function(object, value) standardGeneric("logfc<-"))
setReplaceMethod("logfc", "PhIPData", function(object, value) {
  .replaceAssay(object, "logfc", value)
})

setGeneric("prob<-", function(object, value) standardGeneric("prob<-"))
setReplaceMethod("prob", "PhIPData", function(object, value) {
  .replaceAssay(object, "prob", value)
})

setGeneric("peptideInfo<-", function(object, value) standardGeneric("peptideInfo<-"))
# setReplaceMethod("peptideInfo", "PhIPData", function(object, value){
#   # check # of peptides match with assays
#   if(nrow(value) != nrow(counts(object))){
#     stop("The number of peptides in the annotation differ from the number of peptides in `counts`.")
#   }
#
#   rownames()
#
# })
setGeneric("sampleInfo<-", function(object, value) standardGeneric("sampleInfo<-"))

setGeneric("sampleNames<-", function(object, value) standardGeneric("sampleNames<-"))
setGeneric("peptideNames<-", function(object, value) standardGeneric("peptideNames<-"))


.replaceAssay <- function(object, assay, value){

  # Ensure that target dimensions match
  num_samples <- nrow(sampleInfo(object))
  num_peptides <- nrow(peptideInfo(object))

  if(num_samples == 0 & num_peptides == 0){
    value <- matrix(nrow = num_peptides, ncol = num_samples)
  }
  dim_match <- c(nrow(value) == num_peptides, ncol(value) == num_samples)
  msg <- paste0("Dimensions of ",
                paste0(c("samples", "peptides")[!dim_match], collapse = " and "),
                " do not match.")
  if(sum(dim_match) < 2) { stop(msg) }

  # harmonize names
  rownames(value) <- rownames(peptideInfo(object))
  colnames(value) <- rownames(sampleInfo(object))

  SummarizedExperiment::assays(object)[[assay]] <- DataFrame(value)

  validObject(object)
  object
}

### ==============================================
### Generics
### ==============================================

setMethod("isEmpty", "PhIPData", function(x) {
  isEmpty(counts(x)) & isEmpty(logfc(x)) &
    isEmpty(prob(x)) & isEmpty(peptideInfo(x)) &
    isEmpty(sampleInfo(x))
})

### ==============================================
### Coercion methods
### ==============================================
