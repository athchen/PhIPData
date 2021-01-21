#' @import methods SummarizedExperiment GenomicRanges IRanges edgeR cli
NULL

### PhIPData class  ==============================================
#' PhIPData - A class for PhIP-Seq experiment data
#'
#' @description The \code{PhIPData} class is used to manage the results from
#' phage-immunoprecipitation (PhIP-Seq) experiments. New \code{PhIPData} objects
#' can be created using the homonymous constructor.
#'
#' @return A \code{PhIPData} object
.PhIPData <- setClass("PhIPData", contains = "RangedSummarizedExperiment")


### PhIPData constructor =============================================

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
#' @param metadata a \code{list} object containing experiment-specific metadata
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
                     metadata = list(),
                     .defaultNames = "info", ...) {

  ## Variables defined for convenience
  assays <- c("counts", "logfc", "prob")
  assays_missing <- assays[c(missing(counts), missing(logfc), missing(prob))]

  assay_list <- list(counts = counts, logfc = logfc, prob = prob)
  .defaultNames <- if(length(.defaultNames) == 1) { rep(.defaultNames, 2) } else { .defaultNames }

  ## Check that input dimensions are matched.
  dims <- .checkDims(counts, logfc, prob, peptideInfo, sampleInfo)
  if(is.character(dims)) { stop(dims) }

  ## Get peptide names
  ## if no peptide names are given but there are peptides, assign names of "pep_rownumber"
  peptide_names <- .getPeptideNames(counts, logfc, prob, peptideInfo, .defaultNames[1])
  peptide_names <- if (is.null(peptide_names) & dims[1] != 0) {
    paste0("pep_", 1:dims[1])
    } else { peptide_names }

  ## Get sample names
  ## if no sample names are given but there are samples, assign names of "sample_colnumber"
  sample_names <- .getSampleNames(counts, logfc, prob, sampleInfo, .defaultNames[2])
  sample_names <- if(is.null(sample_names) & dims[2] != 0) { paste0("sample_", 1:dims[2])} else {sample_names}

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
                  prob = assay_list[["prob"]]),
    rowRanges = row_info,
    colData = sample_meta,
    metadata = metadata)

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
    cli::cli_alert_warning(paste0(peptide_warning, "'peptideInfo'."))
    peptide_names <- rownames(peptideInfo)
  } else if(default == "counts") {
    cli::cli_alert_warning(paste0(peptide_warning, "'counts'."))
    peptide_names <- rownames(counts)
  } else if (default== "logfc") {
    cli::cli_alert_warning(paste0(peptide_warning, "'logfc'."))
    peptide_names <- rownames(logfc)
  } else if (default == "prob") {
    cli::cli_alert_warning(paste0(peptide_warning, "'prob'."))
    peptide_names <- rownames(prob)
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
    cli::cli_alert_warning(paste0(sample_warning, "'sampleInfo'."))
    sample_names <- rownames(sampleInfo)
  } else if(default == "counts") {
    cli::cli_alert_warning(paste0(sample_warning, "'counts'."))
    sample_names <- colnames(counts)
  } else if (default == "logfc") {
    cli::cli_alert_warning(paste0(sample_warning, "'logfc'"))
    sample_names <- colnames(logfc)
  } else if (default == "prob") {
    cli::cli_alert_warning(paste0(sample_warning, "'prob'."))
    sample_names <- colnames(prob)
  } else { stop(default_error) }

  sample_names
}

.tidyAssays <- function(assay_list, assays_missing, peptide_names, sample_names){

  assays <- c("counts", "logfc", "prob")
  assays_present <- assays[!assays %in% assays_missing]
  num_missing <- length(assays_missing)

  empty_mat <- matrix(nrow = length(peptide_names),
                      ncol = length(sample_names))
  rownames(empty_mat) <- peptide_names
  colnames(empty_mat) <- sample_names

  for(assay in assays_missing){
    assay_list[[assay]] <- S4Vectors::DataFrame(empty_mat)
  }

  for(assay in assays_present){
    assay_list[[assay]] <- S4Vectors::DataFrame(assay_list[[assay]])
    if(!is.null(peptide_names)) { rownames(assay_list[[assay]]) <- peptide_names }
    if(!is.null(sample_names)) { colnames(assay_list[[assay]]) <- sample_names }
  }

  assay_list
}

.tidyPeptideInfo <- function(peptideInfo, peptide_names){
  peptideInfo <- DataFrame(peptideInfo)

  if(all(c(dim(peptideInfo), length(peptide_names)) == 0)){
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
    if(length(msg) > 0){ cli::cli_alert_warning(msg) }

    if(!all(c("pos_start", "pos_end") %in% colnames(peptideInfo))){
      pep_start <- rep(0, length(peptide_names))
      pep_end <- rep(0, length(peptide_names))
    } else {
      pep_start <- if("pos_start" %in% colnames(peptideInfo)){
        replace(peptideInfo[["pos_start"]], is.na(peptideInfo[, "pos_start"]), 0)
      } else {
        rep(0, length(peptide_names))
      }

      pep_end <- if("pos_end" %in% colnames(peptideInfo)){
        replace(peptideInfo[["pos_end"]], is.na(peptideInfo[, "pos_end"]), 0)
      } else {
        rep(0, length(length(peptide_names)))
      }
    }
  }

  list(pep_meta = pep_meta,
       pep_start = pep_start,
       pep_end = pep_end)
}

### Validity ==============================================
## 1. `counts`, `logfc`, and `prob` must be matrices in a valid PhIPData object.
.checkAssays <- function(x){

  assay_list <- c("counts", "logfc", "prob")

  if(!all(assay_list %in% names(assays(x)))){

    missing_assays <- assay_list[!assay_list %in% names(assays(x))]
    paste0("`counts`, `logfc`, and `prob` assays must be included ",
           "in a PhIPData object. The following assays are missing: ",
           paste(missing_assays, collapse = ", "), ".")

  } else NULL
}

## 2. counts cannot have negative entries
.checkCounts <- function(x){
  counts_df <- as.data.frame(counts(x))
  if(!all(counts_df >= 0 | is.na(counts_df))) {
    "'counts' cannot have negative entries."
  } else NULL

}
## 3. Sample and peptide dimensions must be the same within the object.
.checkObjectDim <- function(x){
  dim_check <- .checkDims(counts(x), logfc(x), prob(x), peptideInfo(x), sampleInfo(x))
  if(is.character(dim_check)) { dim_check } else NULL
}

# 4. sample and peptide names must be identical across all assays and annotation information.
.checkNames <- function(x){

  sample_names <- list(colnames(counts(x)), colnames(logfc(x)),
                    colnames(prob(x)), rownames(sampleInfo(x)))
  peptide_names <- list(rownames(counts(x)), rownames(logfc(x)),
                     rownames(prob(x)), names(peptideInfo(x)))

  match <- c(length(unique(sample_names)) == 1, length(unique(peptide_names)) == 1)

  error <- paste0("Names do not match across ",
                  paste0(c("samples", "peptides")[!match], collapse = " and "), ".")

  if(sum(match) < 2) { error } else { NULL }

}

.validPhIPData <- function(x){
  if(!isEmpty(x)){
    c(.checkAssays(x), .checkCounts(x), .checkObjectDim(x), .checkNames(x))
  }
}

S4Vectors::setValidity2("PhIPData", .validPhIPData)

### ==============================================
### Getters
### ==============================================

#' @export
setGeneric("counts", function(x, ...) standardGeneric("counts"))
setMethod("counts", "PhIPData", function(x)
  SummarizedExperiment::assays(x)[["counts"]])

#' @export
setGeneric("logfc", function(x, ...) standardGeneric("logfc"))
setMethod("logfc", "PhIPData", function(x)
  SummarizedExperiment::assays(x)[["logfc"]])

#' @export
setGeneric("prob", function(x, ...) standardGeneric("prob"))
setMethod("prob", "PhIPData", function(x)
  SummarizedExperiment::assays(x)[["prob"]])

#' @export
setGeneric("peptideInfo", function(x, ...) standardGeneric("peptideInfo"))
setMethod("peptideInfo", "PhIPData", function(x)
  SummarizedExperiment::rowRanges(x))

#' @export
setGeneric("sampleInfo", function(x, ...) standardGeneric("sampleInfo"))
setMethod("sampleInfo", "PhIPData", function(x)
  SummarizedExperiment::colData(x))

### ==============================================
### Setters
### ==============================================

# This `assays` and `assay` replacement functions differs from the
# SummarizedExperiment assays functions in that mismatched names returns
# a valid object rather than an error.
setReplaceMethod("assays", c("PhIPData", "list"), function(x, ..., value) {

  pep_names <- rownames(x)
  sample_names <- colnames(x)

  value <- if(length(value) == 1){
    rownames(value) <- pep_names
    colnames(value) <- sample_names

    value
  } else if(length(value) > 1) {
    lapply(value, function(assay){
      rownames(assay) <- pep_names
      colnames(assay) <- sample_names

      assay
    })
  } else value

  new_object <- callNextMethod()

  # Ensure that counts, logfc, and probs are in the final object
  error <- .checkAssays(new_object)
  if(length(error)){
    stop(error)
  }

  new_object
})


setReplaceMethod("assays", c("PhIPData", "SimpleList"), function(x, ..., value) {

  pep_names <- rownames(x)
  sample_names <- colnames(x)

  value <- if(length(value) == 1){
    rownames(x) <- pep_names
    colnames(x) <- sample_names

    value
  } else if(length(value) > 1) {
    lapply(value, function(assay){
      rownames(assay) <- pep_names
      colnames(assay) <- sample_names

      assay
    })
  } else value

  new_object <- callNextMethod()

  # Ensure that counts, logfc, and probs are in the final object
  error <- .checkAssays(new_object)
  if(length(error)){
    stop(error)
  }

  new_object
})

setReplaceMethod("assay", c("PhIPData", "missing"), function(x, i, ..., value) {

  new_object <- if(!is.null(value)){
    rownames(value) <- rownames(x)
    colnames(value) <- colnames(x)

    callNextMethod()

  } else {
    remaining_assays <- setdiff(assayNames(x), "counts")
    assays(x) <- assays(x)[remaining_assays]

    x
  }

  # Ensure that counts, logfc, and probs are in the final object
  error <- .checkAssays(new_object)
  if(length(error)){
    stop(error)
  }

  new_object
})

setReplaceMethod("assay", c("PhIPData", "numeric"), function(x, i, ..., value) {

  new_object <- if(!is.null(value)){
    rownames(value) <- rownames(x)
    colnames(value) <- colnames(x)

    callNextMethod()

  } else {
    remaining_assays <- assayNames(x)[-i]
    assays(x) <- assays(x)[remaining_assays]

    x
  }

  # Ensure that counts, logfc, and probs are in the final object
  error <- .checkAssays(new_object)
  if(length(error)){
    stop(error)
  }

  new_object
})

setReplaceMethod("assay", c("PhIPData", "character"), function(x, i, ..., value) {

  new_object <- if(!is.null(value)){
    rownames(value) <- rownames(x)
    colnames(value) <- colnames(x)

    callNextMethod()

  } else {
    remaining_assays <- setdiff(assayNames(x), i)
    assays(x) <- assays(x)[remaining_assays]

    x
  }

  # Ensure that counts, logfc, and probs are in the final object
  error <- .checkAssays(new_object)
  if(length(error)){
    stop(error)
  }

  new_object
})

# Convenience functions for standard Assays
#' @export
setGeneric("counts<-", function(x, value) standardGeneric("counts<-"))
setReplaceMethod("counts", "PhIPData", function(x, value) {
  assay(x, "counts") <- value
  x
})

#' @export
setGeneric("logfc<-", function(x, value) standardGeneric("logfc<-"))
setReplaceMethod("logfc", "PhIPData", function(x, value) {
  assay(x, "logfc") <- value
  x
})

#' @export
setGeneric("prob<-", function(x, value) standardGeneric("prob<-"))
setReplaceMethod("prob", "PhIPData", function(x, value) {
  assay(x, "prob") <- value
  x
})

#' @export
setGeneric("peptideInfo<-", function(x, value) standardGeneric("peptideInfo<-"))
setReplaceMethod("peptideInfo", "PhIPData", function(x, value){
  # check # of peptides match with assays
  if(nrow(value) != nrow(counts(x))){
    stop("The number of peptides in the annotation differ from the number of peptides in `counts`.")
  }

  # get peptide names and sample names from existing x
  rownames(value) <- peptide_names <- dimnames(x)[[1]]

  new_pepInfo <- .tidyPeptideInfo(value, peptide_names)

  pep_meta <- new_pepInfo[["pep_meta"]]
  pep_start <- new_pepInfo[["pep_start"]]
  pep_end <- new_pepInfo[["pep_end"]]

  row_info <- GenomicRanges::GRanges(seqnames = peptide_names,
                                     ranges = IRanges::IRanges(start = pep_start, end = pep_end))
  if(!S4Vectors::isEmpty(pep_meta)){ mcols(row_info) <- pep_meta }

  # Why does setting a new rowRanges erase all my row names for the x?!?!?!?!?
  rowRanges(x) <- row_info
  rownames(x) <- peptide_names

  validObject(x)

  x
})

#' @export
setGeneric("sampleInfo<-", function(x, value) standardGeneric("sampleInfo<-"))
setReplaceMethod("sampleInfo", "PhIPData", function(x, value){
  # check # of samples match with assays
  if(nrow(value) != ncol(counts(x))){
    stop("The number of samples in the annotation differ from the number of samples in `counts`.")
  }

  # get peptide names and sample names from existing object
  rownames(value) <- sample_names <- dimnames(x)[[2]]

  sample_meta <- S4Vectors::DataFrame(value, row.names = sample_names)
  colData(x) <- sample_meta

  validObject(x)

  x
})

### ==============================================
### Generics
### ==============================================

# setMethod("show", "PhIPData", function(object) {
#
# })

setMethod("isEmpty", "PhIPData", function(x) {
  isEmpty(counts(x)) & isEmpty(logfc(x)) &
    isEmpty(prob(x)) & isEmpty(peptideInfo(x)) &
    isEmpty(sampleInfo(x))
})

### ==============================================
### Coercion methods
### ==============================================

# DGEList to PhIPData
setAs("DGEList", "PhIPData", function(from){

  count_mat <- if(is.null(from[["counts"]])) { S4Vectors::DataFrame() } else { from[["counts"]] }
  sample_info <- if(is.null(from[["samples"]])) { S4Vectors::DataFrame() } else { from[["samples"]] }
  pep_info <- if(is.null(from[["genes"]])) { S4Vectors::DataFrame() } else { from[["genes"]] }

  PhIPData(counts = count_mat, peptideInfo = pep_info, sampleInfo = sample_info)
})

# PhIPData to DGEList
setAs("PhIPData", "DGEList", function(from){
  count_mat <- counts(from)
  sample_info <- sampleInfo(from)
  peptide_info <- peptideInfo(from)
  group_labs <- if("group" %in% colnames(sample_info)){ sample_info$group } else NULL

  edgeR::DGEList(counts = count_mat,
          samples = sample_info[colnames(sample_info) != "group"],
          group = group_labs,
          genes = peptide_info)
})

# PhIPData to list
setAs("PhIPData", "list", function(from){

  list(assays = assays(from),
       peptideInfo = peptideInfo(from),
       sampleInfo = sampleInfo(from),
       metadata = metadata(from))
})

setAs("list", "PhIPData", function(from){

  PhIPData(counts = from$assays[["counts"]],
           logfc = from$assays[["logfc"]],
           prob = from$assays[["prob"]],
           peptideInfo = from$peptideInfo,
           sampleInfo = from$sampleInfo,
           metadata = from$metadata)
})


# PhIPData to List
setAs("PhIPData", "List", function(from){

  List(assays = assays(from),
       peptideInfo = peptideInfo(from),
       sampleInfo = sampleInfo(from),
       metadata = metadata(from))
})

setAs("List", "PhIPData", function(from){

  PhIPData(counts = from[["assays"]][["counts"]],
           logfc = from[["assays"]][["logfc"]],
           prob = from[["assays"]][["prob"]],
           peptideInfo = from[["peptideInfo"]],
           sampleInfo = from[["sampleInfo"]],
           metadata = from[["metadata"]])
})
