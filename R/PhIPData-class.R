#' @import methods SummarizedExperiment GenomicRanges IRanges
NULL

### ==============================================
### PhIPData class
### ==============================================

#' PhIPData - A class for PhIP-Seq experiment data
#'
#' @description The \code{PhIPData} class is used to manage the results from
#' phage-immunoprecipitation (PhIP-Seq) experiments. New \code{PhIPData} objects
#' can be created using the helper function \code{makePhIPData} (recommended) or the homonymous constructor.
#'
#' @return A \code{PhIPData} object
.PhIPData <- setClass("PhIPData", contains = "SummarizedExperiment")

### ==============================================
### PhIPData constructor
### ==============================================

#' Constructor for PhIPData object.
#' @export
#' @param counts matrix of integer read counts
#' @param logfc matrix of log10 fold changes
#' @param prob matrix of p-values or posterior probabilities for estimated enrichment
#' @param peptideInfo data frame of peptide information, must contain columns `pep_id`, `pos_start`, and `pos_end`
#' @param sampleInfo data frame of sample information.
PhIPData <- function(counts, logfc, prob, peptideInfo, sampleInfo, ...) {
  # Extract metadata information that is not position-related.
  pep_meta <- peptideInfo[,!colnames(peptideInfo) %in% c("pep_id", "pos_start", "pos_end")]

  row_info <- GenomicRanges::GRanges(seqnames = peptideInfo$pep_id,
                      ranges = IRanges::IRanges(start = peptideInfo$pos_start,
                                       end = peptideInfo$pos_end))
  mcols(row_info) <- pep_meta

  ## Make RangedSummarizedExperiment
  se_object <- SummarizedExperiment::SummarizedExperiment(list(counts = counts, logfc = logfc, prob = prob),
                         rowRanges = row_info,
                         colData = sampleInfo)

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
