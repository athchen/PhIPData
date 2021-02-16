#' @import SummarizedExperiment methods S4Vectors
#' @include defineBeads.R
NULL

### PhIPData class  ==============================================
#' The PhIPData class
#'
#' @description
#' The \code{PhIPData} class is a matrix-like container designed to organize
#' results from phage-immunoprecipitation (PhIP-Seq) experiments. Rows in
#' PhIPData objects represent peptides and columns represent samples. Each
#' object contains at least three assays:
#'
#' \itemize{
#'     \item{\code{counts}:} a matrix of raw read counts,
#'     \item{\code{logfc}:} a matrix of log10 estimated fold-change in
#'          comparison to beads-only samples,
#'     \item{\code{prob}:} a matrix of probabilities associated with whether
#'          a sample has an enriched antibody response for a peptide.
#' }
#'
#' The \code{PhIPData} class extends the
#' \linkS4class{RangedSummarizedExperiment} class, so methods documented
#' in \linkS4class{RangedSummarizedExperiment} and
#' \linkS4class{SummarizedExperiment} also work on \code{PhIPData} objects.
#'
#' @details
#' Rows of \code{PhIPData} objects correspond to peptides of interest and are
#' organized in \link{GRanges} or \link{GRangesList} objects. Though originally
#' designed for genomic ranges, the sequence name and genomic range information
#' in \link{GRanges} objects can be replaced with peptide names and amino acid
#' positions, respectively. If no peptide names are given, peptides are given
#' the names of \code{pep_rownum}. Peptide positions are specified by columns
#' \code{pos_start} and \code{pos_end} in the \code{peptideInfo} argument of the
#' constuctor. Missing position information is set to 0. Additional peptide
#' annotation can also be stored in \link{GRanges} objects and can be used
#' to subset \code{PhIPData} objects as shown below.
#'
#' Columns of \code{PhIPData} objects represent samples. Sample metadata
#' are stored in a \link{DataFrame} and can be accessed as shown below. If no
#' sample names are specified, samples are given default names of
#' \code{sample_colnum}.
#'
#' Unlike \link{RangedSummarizedExperiment}/\link{SummarizedExperiment} objects,
#' \code{PhIPData} objects must contain \code{counts}, \code{logfc}, \code{prob}.
#' If any of the three assays are missing when the constructor is called, an
#' empty matrix of the same names and dimensions is initialized for that assay.
#' Sample and peptide names are harmonized across assays and annotation during
#' construction and replacement.
#'
#' Though `counts` typically contain integer values for the number of reads
#' aligned to each peptide, `PhIPData` only requires that stored values are
#' non-negative numeric values. Pseudocounts or non-integer count values can also
#' be stored in the `counts` assay.
#'
#' @section Constructor:
#' \code{PhIPData} objects are constructed using the homonymous function and
#' arguments as described above. Any \code{PhIPData} object can be created
#' so long as peptide and sample identifiers (or lack thereof) are specified
#' via any of the parameters.
#'
#' @seealso
#'      \code{\link{PhIPData-methods}} for accessors and modifiers for PhIPData
#'      components.
#'      \linkS4class{SummarizedExperiment}
.PhIPData <- setClass("PhIPData", contains = "RangedSummarizedExperiment")

### PhIPData constructor =============================================
#' @rdname PhIPData-class
#'
#' @param counts a \code{matrix}, \code{data.frame}, or
#'     \code{\linkS4class{DataFrame}}of \strong{integer} read counts.
#' @param logfc a \code{matrix}, \code{data.frame}, or
#'     \code{\linkS4class{DataFrame}} of log10 estimated fold changes.
#' @param prob a \code{matrix}, \code{data.frame}, or
#'     \code{\linkS4class{DataFrame}} of probability values (p-values or
#'     posterior probabilities) for enrichment estimates.
#' @param peptideInfo a \code{data.frame} or \code{\linkS4class{DataFrame}} of
#'    peptide information.
#' @param sampleInfo a \code{data.frame} or \code{\linkS4class{DataFrame}} of
#'     additional sample information.
#' @param metadata a \code{list} object containing experiment-specific metadata.
#' @param .defaultNames vector of names to use when sample and peptide
#'     identifiers disagree across the metadata and the \code{counts},
#'     \code{logfc}, and \code{prob} matrices. If \code{.defaultNames} is of
#'     length 1, the same source is used for both peptide and sample
#'     identifiers. If \code{.defaultNames} is longer than 2, the first and
#'     second elements correspond to the names for peptides and samples,
#'     respectively.
#'
#'     Valid options are:
#'     \itemize{
#'         \item{"info": }{names should be taken from the \code{SampleInfo} or
#'             \code{peptideInfo} objects.}
#'         \item{"counts": }{names should be taken from the row/column names
#'             of the \code{counts} object.}
#'         \item{"logfc": }{names should be taken from the row/column names of
#'             the \code{logfc} object.}
#'         \item{"prob": }{names should be taken from the row/column names of
#'             the \code{prob} object.}
#'     }
#'
#' @return A \code{PhIPData} object.
#'
#' @examples
#' ## Construct a new PhIPData object
#' counts_dat <- matrix(sample(1:1e6, 25, replace = TRUE), nrow = 5)
#' logfc_dat <- matrix(rnorm(25, 0, 10), nrow = 5)
#' prob_dat <- matrix(rbeta(25, 1, 1), nrow = 5)
#'
#' peptide_meta <- data.frame(pos_start = 1:5,
#'       pos_end = 6:10,
#'       species = c(rep("HIV", 3), rep("EBV", 2)))
#' sample_meta <- data.frame(gender = sample(c("M", "F"), 5, TRUE),
#'      group = sample(c("ctrl", "trt", "beads"), 5, TRUE))
#' exp_meta<- list(date_run = as.Date("2021/01/20"),
#'       reads_per_sample = colSums(counts_dat))
#'
#' rownames(counts_dat) <- rownames(logfc_dat) <-
#'      rownames(prob_dat) <- rownames(peptide_meta) <-
#'      paste0("pep_", 1:5)
#' colnames(counts_dat) <- colnames(logfc_dat) <-
#'      colnames(prob_dat) <- rownames(sample_meta) <-
#'      paste0("sample_", 1:5)
#'
#' phip_obj <- PhIPData(counts_dat, logfc_dat, prob_dat,
#'      peptide_meta, sample_meta, exp_meta)
#' phip_obj
#'
#' @export
#' @importClassesFrom edgeR DGEList
#' @importFrom S4Vectors DataFrame isEmpty setValidity2
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom cli cli_alert_warning
PhIPData <- function(counts = S4Vectors::DataFrame(),
                     logfc = S4Vectors::DataFrame(),
                     prob = S4Vectors::DataFrame(),
                     peptideInfo = S4Vectors::DataFrame(),
                     sampleInfo = S4Vectors::DataFrame(),
                     metadata = list(),
                     .defaultNames = "info") {

  ## Variables defined for convenience
  assays <- c("counts", "logfc", "prob")
  assays_missing <- assays[c(missing(counts), missing(logfc), missing(prob))]

  assay_list <- list(counts = counts, logfc = logfc, prob = prob)
  .defaultNames <- if(length(.defaultNames) == 1) {
    rep(.defaultNames, 2)
  } else {
    .defaultNames
  }

  ## Check that input dimensions are matched.
  dims <- .checkDims(counts, logfc, prob, peptideInfo, sampleInfo)
  if(is.character(dims)) { stop(dims) }

  ## Get peptide names. If no peptide names are given but there are peptides,
  ## assign names of "pep_rownumber"
  peptide_names <- .getPeptideNames(counts, logfc, prob, peptideInfo, .defaultNames[1])
  peptide_names <- if(is.null(peptide_names) & dims[1] != 0) {
    paste0("pep_", seq_len(dims[1]))
    } else { peptide_names }

  ## Get sample names. If no sample names are given but there are samples,
  ## assign names of "sample_colnumber"
  sample_names <- .getSampleNames(counts, logfc, prob, sampleInfo, .defaultNames[2])
  sample_names <- if(is.null(sample_names) & dims[2] != 0) {
    paste0("sample_", seq_len(dims[2]))
  } else { sample_names }

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
  sample_meta <- if(missing(sampleInfo)){
    S4Vectors::DataFrame(group = rep(NA, length = length(sample_names)),
                         row.names = sample_names)
  } else {

    ## Add 'group' column if it is not present
    if(!"group" %in% colnames(sampleInfo)){
      cli_alert_warning("No 'group' column in sampleInfo. Adding empty column.")
      sampleInfo$group <- NA
    }

    ## Warn if there are no beads-only samples
    if(!getBeadsName() %in% sampleInfo$group){
      cli_alert_info("No beads-only samples present in the PhIPData object.")
    }

    S4Vectors::DataFrame(sampleInfo, row.names = sample_names)
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

  fixed_row <- if(length(unique(row_dims)) == 1) {
    unique(row_dims)
  } else { unique(row_dims[row_dims != 0]) }
  fixed_col <- if(length(unique(col_dims)) == 1) {
    unique(col_dims)
  } else { unique(col_dims[col_dims != 0]) }

  match <- c(length(fixed_row) == 1, length(fixed_col) == 1)
  error <- paste0("The number of ",
                  paste0(c("samples", "peptides")[! match], collapse = " and "),
                  " differs across inputs.")

  if(sum(match) < 2){ error } else { c(fixed_row, fixed_col) }
}

.getPeptideNames <- function(counts, logfc, prob, peptideInfo, default){

  peptide_names <- list(rownames(counts), rownames(logfc),
                        rownames(prob), rownames(peptideInfo))
  non_missing_len <- vapply(peptide_names, length, numeric(1))
  peptide_names <- unique(peptide_names[non_missing_len != 0])

  peptide_warning <- paste0("Peptide names are not identical across inputs. ",
                            "Using peptide names from ")
  default_error <- paste0("Invalid '.defaultNames' supplied. Valid ",
                          "'.defaultNames' options are 'info', 'counts', ",
                          "'logfc', or 'prob'.")

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

  sample_names <- list(colnames(counts), colnames(logfc),
                       colnames(prob), rownames(sampleInfo))
  non_missing_len <- vapply(sample_names, length, numeric(1))
  sample_names <- unique(sample_names[non_missing_len != 0])

  sample_warning <- paste0("Sample names are not identical across inputs. ",
                           "Using sample names from ")
  default_error <- paste0("Invalid '.defaultNames' supplied. Valid ",
                          "'.defaultNames' options are: 'info', 'counts', ",
                          "'logfc', or 'prob'.")

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

  ## Create empty matrix missing assays
  empty_mat <- matrix(nrow = length(peptide_names),
                      ncol = length(sample_names))
  rownames(empty_mat) <- peptide_names
  colnames(empty_mat) <- sample_names
  for(assay in assays_missing){
    assay_list[[assay]] <- S4Vectors::DataFrame(empty_mat)
  }

  ## Coerce all assays into DataFrames
  assay_list <- lapply(assay_list, S4Vectors::DataFrame)

  ## Correct names
  for(assay in assays){
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
    pep_meta <- S4Vectors::DataFrame(peptideInfo[, !colnames(peptideInfo) %in%
                                                   c("pos_start", "pos_end"),
                                                 drop = FALSE])

    warnings <- c(no_start = (!"pos_start" %in% colnames(peptideInfo)),
                  no_end = (!"pos_end" %in% colnames(peptideInfo)))
    warnings["missing_start"] <- if(!warnings["no_start"]){
      any(is.na(peptideInfo[, "pos_start"])) } else { TRUE }
    warnings["missing_end"] <- if(!warnings["no_end"]){
      any(is.na(peptideInfo[, "pos_end"]))
    } else { TRUE }

    msg <- if (sum(warnings[c("no_start", "no_end")]) > 0) {
      paste0("Missing peptide ",
             paste(c("start", "end")[warnings[c("no_start", "no_end")]],
                   collapse = " and "),
             " position information.")
    } else NULL

    msg <- if(any(warnings)) {
      paste0(msg, " Replacing missing values with 0.")
    } else { msg }
    if(length(msg) > 0){ cli::cli_alert_warning(msg) }

    ## Add pos_start and pos_end columns if not present, set to 0 if missing
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
  error <- character()

  if(all(dim(x) != 0)) {
    counts_pos <- vapply(counts(x), function(col) col >= 0 | is.na(col),
                         logical(nrow(x)))
    if(!all(counts_pos)) {
      msg <- "cannot have negative entries"
      error <- c(error, msg)
    }

    error <- paste0(error, collapse = " and ")

    if(error != ""){
      paste0("'counts' ", error, ".")
    } else NULL
  } else NULL
}

## 3. Sample and peptide dimensions must be the same within the object.
.checkObjectDim <- function(x){
  dim_check <- .checkDims(counts(x), logfc(x), prob(x), peptideInfo(x), sampleInfo(x))
  if(is.character(dim_check)) { dim_check } else NULL
}

## 4. sample and peptide names must be identical across all assays and annotation information.
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

.checkSampleInfo <- function(x){
  if(!"group" %in% colnames(sampleInfo(x))){
    "'group' column must be present in sample information."
  } else NULL
}

.validPhIPData <- function(x){
  if(!isEmpty(x)){
    c(.checkAssays(x), .checkCounts(x), .checkObjectDim(x), .checkNames(x),
      .checkSampleInfo(x))
  }
}

S4Vectors::setValidity2("PhIPData", .validPhIPData)

### Getters ==============================================
#' @name PhIPData-methods
#' @title Accessing and Modifying Information in PhIPData objects
#'
#' @description Methods to extract and modify \code{assay(s)}(including
#' convenient functions for \code{counts}, \code{logfc}, and \code{prob}),
#' \code{sampleInfo}, \code{peptideInfo}, and \code{metadata}.
#'
#' @details In addition to the functions detailed in
#' \linkS4class{RangedSummarizedExperiment}, the \code{PhIPData} class includes
#' conveniently named functions to quickly access and modify frequently used
#' components of PhIPData objects.
#'
#' Replacement functions ensure that names of the replacement object are matched
#' with the names of the \code{PhIPData} object.
#'
#' Since packages for identifying differential expression in RNA-seq experiments
#' are frequently used for estimating fold-changes for peptide enrichments,
#' the class also includes coercion methods to and from \linkS4class{DGEList}s.
#'
#' @section Available methods:
#' In the following code snippets, \code{x} is a \linkS4class{PhIPData} object,
#' \code{value} is a matrix-like object with the same dimensions as \code{x},
#' and \code{...} are further arguments passed to \code{\link{assay}} (for the getter) or \code{\link{assay<-}} (for the setter).
#' \describe{
#' \item{\code{counts(x, ...)}, \code{counts(x, ...) <- value}:}{
#' Get or set a matrix of raw read counts
#' }
#' \item{\code{logfc(x, ...)}, \code{logfc(x, ...) <- value}:}{
#' Get or set a matrix of log10 estimated fold changes (in comparison to beads-only samples)
#' }
#' \item{\code{prob(x, ...)}, \code{pob(x, ...) <- value}:}{
#' Get or set a matrix of probabilities associated with whether
#'          a sample has an enriched antibody response for a peptide.
#' }
#' }
#'
#' @examples
#' example("PhIPData")
#'
#' replacement_dat <- matrix(1L, nrow = 5, ncol = 5)
#'
#' ## SummarizedExperiment Accessors and Setters
#' assays(phip_obj)
#' assays(phip_obj)$counts <- replacement_dat
#' assay(phip_obj, "logfc")
#' assay(phip_obj, "logfc") <- replacement_dat
#'
#' ## counts
#' counts(phip_obj)
#' counts(phip_obj) <- counts_dat
#'
#' ## logfc
#' logfc(phip_obj)
#' logfc(phip_obj) <- logfc_dat
#'
#' ## prob
#' prob(phip_obj)
#' prob(phip_obj) <- replacement_dat
#'
#' ## coercion functions
#' as(phip_obj, "DGEList")
#' as(phip_obj, "List")
#' as(phip_obj, "list")
#'
#' @param object A \code{PhIPData} object
#' @param x A \code{PhIPData} object
#' @param i A \code{numeric}, {character}
#' @param withDimnames Parameter for
#'      \linkS4class{RangedSummarizedExperiment} class functions. Overrided
#'      since row/column names are automatically synced within each object.
#' @param ... parameters for \code{\link[SummarizedExperiment]{assays}},
#'      which are typically not needed.
#' @param value A \code{matrix}, \code{data.frame}, or \linkS4class{DataFrame}
#'      of the same dimensions (not necessarily the same names)
#'
#' @return Accessors: a \link{DataFrame} object
#' @return Setters: a \code{PhIPData} object
#'
#' @seealso
#' \code{\link[SummarizedExperiment]{assays}} for
#' \linkS4class{SummarizedExperiment} operations.
NULL

#' @export
#' @rdname PhIPData-methods
#' @importFrom BiocGenerics counts
setMethod("counts", "PhIPData", function(object, ...)
  SummarizedExperiment::assay(object, "counts", ...))

#' @export
#' @rdname PhIPData-methods
setGeneric("logfc", function(object, ...) standardGeneric("logfc"))

#' @export
#' @rdname PhIPData-methods
setMethod("logfc", "PhIPData", function(object, ...)
  SummarizedExperiment::assay(object, "logfc", ...))


#' @export
#' @rdname PhIPData-methods
setGeneric("prob", function(object, ...) standardGeneric("prob"))

#' @export
#' @rdname PhIPData-methods
setMethod("prob", "PhIPData", function(object, ...)
  SummarizedExperiment::assay(object, "prob", ...))

#' @export
#' @rdname PhIPData-methods
setGeneric("peptideInfo", function(object, ...) standardGeneric("peptideInfo"))

#' @export
#' @rdname PhIPData-methods
setMethod("peptideInfo", "PhIPData", function(object, ...)
  SummarizedExperiment::rowRanges(object, ...))


#' @export
#' @rdname PhIPData-methods
setGeneric("sampleInfo", function(object, ...) standardGeneric("sampleInfo"))

#' @export
#' @rdname PhIPData-methods
setMethod("sampleInfo", "PhIPData", function(object, ...)
  SummarizedExperiment::colData(object, ...))

### Setters ==============================================
# This `assays` and `assay` replacement functions differs from the
# SummarizedExperiment assays functions in that mismatched names returns
# a valid object rather than an error.
#' @export
#' @rdname PhIPData-methods
setReplaceMethod("assays", c("PhIPData", "list"), function(x, ..., value) {

  pep_names <- rownames(x)
  sample_names <- colnames(x)

  value <- if(length(value) == 1){
    rownames(value) <- pep_names
    colnames(value) <- sample_names

    S4Vectors::DataFrame(value)
  } else if(length(value) > 1) {
    lapply(value, function(assay){
      rownames(assay) <- pep_names
      colnames(assay) <- sample_names

      S4Vectors::DataFrame(assay)
    })
  } else value

  new_object <- callNextMethod()

  # Ensure that counts, logfc, and probs are in the final object
  error <- .checkAssays(new_object)
  if(length(error)){
    stop(error)
  }

  validObject(new_object)
  new_object
})

#' @export
#' @rdname PhIPData-methods
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

  validObject(new_object)
  new_object
})

#' @export
#' @rdname PhIPData-methods
setReplaceMethod("assay", c("PhIPData", "missing"),
                 function(x, i, withDimnames = TRUE, ..., value) {

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

  validObject(new_object)
  new_object
})

#' @export
#' @rdname PhIPData-methods
setReplaceMethod("assay", c("PhIPData", "numeric"),
                 function(x, i, withDimnames = TRUE, ..., value) {

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

  validObject(new_object)
  new_object
})

#' @export
#' @rdname PhIPData-methods
setReplaceMethod("assay", c("PhIPData", "character"),
                 function(x, i, withDimnames = TRUE, ..., value) {

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

  validObject(new_object)
  new_object
})

# Convenience functions for standard Assays
#' @export
#' @rdname PhIPData-methods
#' @importFrom BiocGenerics "counts<-"
setReplaceMethod("counts", c("PhIPData", "ANY"), function(object, ..., value) {
  assay(object, "counts") <- value
  object
})

#' @export
#' @rdname PhIPData-methods
setGeneric("logfc<-", function(object, ..., value) standardGeneric("logfc<-"))

#' @export
#' @rdname PhIPData-methods
setReplaceMethod("logfc", "PhIPData", function(object, ..., value) {
  assay(object, "logfc") <- value
  object
})


#' @export
#' @rdname PhIPData-methods
setGeneric("prob<-", function(object, ..., value) standardGeneric("prob<-"))

#' @export
#' @rdname PhIPData-methods
setReplaceMethod("prob", "PhIPData", function(object, ..., value) {
  assay(object, "prob") <- value
  object
})

#' @export
#' @rdname PhIPData-methods
setGeneric("peptideInfo<-", function(object, value) standardGeneric("peptideInfo<-"))

#' @export
#' @rdname PhIPData-methods
setReplaceMethod("peptideInfo", "PhIPData", function(object, value){
  # check # of peptides match with assays
  if(nrow(value) != nrow(counts(object))){
    stop("The number of peptides in the annotation differ from the number of peptides in `counts`.")
  }

  # get peptide names and sample names from existing x
  rownames(value) <- peptide_names <- dimnames(object)[[1]]

  new_pepInfo <- .tidyPeptideInfo(value, peptide_names)

  pep_meta <- new_pepInfo[["pep_meta"]]
  pep_start <- new_pepInfo[["pep_start"]]
  pep_end <- new_pepInfo[["pep_end"]]

  row_info <- GenomicRanges::GRanges(seqnames = peptide_names,
                                     ranges = IRanges::IRanges(start = pep_start, end = pep_end))
  if(!S4Vectors::isEmpty(pep_meta)){ mcols(row_info) <- pep_meta }

  # Why does setting a new rowRanges erase all my row names for the x?!?!?!?!?
  rowRanges(object) <- row_info
  rownames(object) <- peptide_names

  validObject(object)

  object
})

#' @export
#' @rdname PhIPData-methods
setGeneric("sampleInfo<-", function(object, ..., value) standardGeneric("sampleInfo<-"))

#' @export
#' @rdname PhIPData-methods
#' @importFrom cli cli_alert_warning cli_alert_info
setReplaceMethod("sampleInfo", "PhIPData", function(object, value){
  # check # of samples match with assays
  if(nrow(value) != ncol(counts(object))){
    stop("The number of samples in the annotation differ from the number of samples in `counts`.")
  }

  # get peptide names and sample names from existing object
  rownames(value) <- sample_names <- dimnames(object)[[2]]

  sample_meta <- S4Vectors::DataFrame(value, row.names = sample_names)

  ## Add 'group' column if it is not present
  if(!"group" %in% colnames(sample_meta)){
    cli::cli_alert_warning("No 'group' column in sampleInfo. Adding empty column.")
    sample_meta$group <- NA
  }

  ## Warn if there are no beads-only samples
  if(!getBeadsName() %in% sample_meta$group){
    cli::cli_alert_info("No beads-only samples present in the PhIPData object.")
  }

  colData(object) <- sample_meta

  validObject(object)

  object
})

### Subsetters ==============================================

### Generics ==============================================
setMethod("show", "PhIPData", function(object) {
  callNextMethod()
  beads_name <- getBeadsName()
  num_beads <- sum(sampleInfo(object)$group == beads_name, na.rm = TRUE)
  cat("beads-only name(", num_beads, "): ", beads_name,
      sep = "")
})

setMethod("isEmpty", "PhIPData", function(x) {
  isEmpty(counts(x)) & isEmpty(logfc(x)) &
    isEmpty(prob(x)) & isEmpty(peptideInfo(x)) &
    isEmpty(sampleInfo(x))
})

### Coercion methods ==============================================
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
