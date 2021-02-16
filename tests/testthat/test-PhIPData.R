context("Base PhIPData API functions work as expected.")

source("setup.R")

# Test missing params ----------------------------------------
test_that("valid PhIPData objects are created when there are missing inputs.", {

  # No assay, no metadata
  expect_is(PhIPData(), "PhIPData")

  # No assay, only metadata
  expect_is(PhIPData(peptideInfo = virscan_info), "PhIPData")
  expect_is(PhIPData(sampleInfo = sampleInfo), "PhIPData")

  # only one assay/input is present, no metadata
  expect_is(PhIPData(counts = counts), "PhIPData")
  expect_is(PhIPData(logfc = logfc), "PhIPData")
  expect_is(PhIPData(prob = prob), "PhIPData")

  # one assay and metadata
  expect_is(PhIPData(counts = counts, sampleInfo = sampleInfo), "PhIPData")
  expect_is(PhIPData(counts = counts, peptideInfo = virscan_info), "PhIPData")
  expect_is(PhIPData(counts = counts, peptideInfo = virscan_info,
                     sampleInfo = sampleInfo), "PhIPData")

  # only two assays are present, no metadata
  expect_is(PhIPData(counts = counts, logfc = logfc), "PhIPData")
  expect_is(PhIPData(counts = counts, prob = prob), "PhIPData")
  expect_is(PhIPData(logfc = logfc, prob = prob), "PhIPData")

  # only two assays and metadata
  expect_is(PhIPData(counts = counts, logfc = logfc,
                     sampleInfo = sampleInfo), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc,
                     peptideInfo = virscan_info), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc,
                     peptideInfo = virscan_info, sampleInfo = sampleInfo),
            "PhIPData")

  # all 3 assays, no metadata
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob), "PhIPData")

  # all 3 assays and metadata
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     sampleInfo = sampleInfo), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     peptideInfo = virscan_info), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     peptideInfo = virscan_info, sampleInfo = sampleInfo),
            "PhIPData")

  # peptide data is missing positional information
  nopos_pepInfo <- dplyr::select(virscan_info, -contains("pos"))
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     peptideInfo = nopos_pepInfo, sampleInfo = sampleInfo),
            "PhIPData")
})

# Test name-fixing code ----------------------------------------
# Function to quickly compare names
expect_names <- function(object, list){
  names <- dimnames(object)
  return(all(names[[1]] == list[[1]]) & all(names[[2]] == list[[2]]))
}

test_that("dimension names are set correctly when inputs are mismatched", {
  assay_list <- list(counts = counts, logfc = logfc, prob = prob)

  for(assay in c("counts","prob")){
    rownames(assay_list[[assay]]) <-
      paste0(assay, "_", rownames(assay_list[[assay]]))
    colnames(assay_list[[assay]]) <-
      paste0(assay, "_", colnames(assay_list[[assay]]))
  }
  colnames(assay_list[["logfc"]]) <-
    paste0("logfc_", colnames(assay_list[["logfc"]]))

  # Check defaults work
  expect_true(expect_names(PhIPData(assay_list[["counts"]],
                                    assay_list[["logfc"]],
                                    assay_list[["prob"]],
                                    virscan_info,
                                    sampleInfo),
                           list(rownames(peptideInfo),
                                rownames(sampleInfo))))

  # Check single default input works as expected
  expect_true(expect_names(PhIPData(assay_list[["counts"]],
                                    assay_list[["logfc"]],
                                    assay_list[["prob"]],
                                    virscan_info,
                                    sampleInfo,
                                    .defaultNames = "counts"),
                           dimnames(assay_list[["counts"]])))

  expect_true(expect_names(PhIPData(assay_list[["counts"]],
                                    assay_list[["logfc"]],
                                    assay_list[["prob"]],
                                    virscan_info,
                                    sampleInfo,
                                    .defaultNames = "logfc"),
                           dimnames(assay_list[["logfc"]])))

  expect_true(expect_names(PhIPData(assay_list[["counts"]],
                                    assay_list[["logfc"]],
                                    assay_list[["prob"]],
                                    virscan_info,
                                    sampleInfo,
                                    .defaultNames = "prob"),
                           dimnames(assay_list[["prob"]])))

  # check double default and extra default input do not generate errors
  expect_true(expect_names(PhIPData(assay_list[["counts"]],
                                    assay_list[["logfc"]],
                                    assay_list[["prob"]],
                                    virscan_info,
                                    sampleInfo,
                                    .defaultNames = c("counts", "logfc",
                                                      "extraneous")),
                           list(rownames(assay_list[["counts"]]),
                                colnames(assay_list[["logfc"]]))))
})

# Test error-causing inputs ----------------------------------------
test_that("invalid inputs return proper errors", {

  # negative counts
  expect_error(PhIPData(counts = matrix(-15:14, nrow = 5)),
               "'counts' cannot have negative entries.")

  # incompatible dimensions
  expect_error(PhIPData(counts = counts,
                        logfc = DataFrame(matrix(nrow = 5, ncol = 5))),
               "The number of .* differs across inputs")

  # invalid default behavior for mismatched names
  # (only a problem names are mismatched)
  rownames(logfc) <- paste0("logfc_", rownames(logfc))
  expect_error(PhIPData(counts, logfc, prob,
                        virscan_info, sampleInfo,
                        .defaultNames = "test"),
               "Invalid '.defaultNames' supplied.")

  colnames(logfc) <- paste0("logfc_", colnames(logfc))
  expect_error(PhIPData(counts, logfc, prob,
                        virscan_info, sampleInfo,
                        .defaultNames = "test"),
               "Invalid '.defaultNames' supplied.")
})

# Test getters ----------------------------------------
test_that("getter functions return objects of the expected class", {
  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo, peptideInfo = virscan_info)

  expect_is(counts(phip_obj), "DataFrame")
  expect_is(logfc(phip_obj), "DataFrame")
  expect_is(prob(phip_obj), "DataFrame")
  expect_is(peptideInfo(phip_obj), "GRanges")
  expect_is(sampleInfo(phip_obj), "DataFrame")
})

# Test setters ----------------------------------------
# 1. Check that replacement works and
# 2. Replacement fixes mismatched sample/peptide names.
test_that("setter functions change the object as desired", {
  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo, peptideInfo = virscan_info)

  # Check assay replacement; newinfo has different names
  replacement_matrix <- matrix(sample(1:1e6, n_samples*n_peptides,
                                      replace = TRUE), nrow = n_peptides)

  expect_error(assays(phip_obj) <- list(assay_1 = replacement_matrix,
                                        assay_2 = replacement_matrix),
               paste0("`counts`, `logfc`, and `prob` assays must be included ",
                      "in a PhIPData object. The following assays are ",
                      "missing: counts, logfc, prob."))
  assays(phip_obj) <- list(counts = replacement_matrix,
                           logfc = replacement_matrix,
                           prob = replacement_matrix)
  expect_equal(unname(as.matrix(counts(phip_obj))), replacement_matrix)
  expect_equal(unname(as.matrix(logfc(phip_obj))), replacement_matrix)
  expect_equal(unname(as.matrix(prob(phip_obj))), replacement_matrix)

  assay(phip_obj) <- counts
  expect_equal(counts(phip_obj), S4Vectors::DataFrame(counts))
  assay(phip_obj, 2) <- logfc
  expect_equal(logfc(phip_obj), S4Vectors::DataFrame(logfc))
  assay(phip_obj, "prob") <- prob
  expect_equal(prob(phip_obj), S4Vectors::DataFrame(prob))

  counts(phip_obj) <- logfc(phip_obj) <- prob(phip_obj) <- replacement_matrix
  expect_equal(unname(as.matrix(counts(phip_obj))), replacement_matrix)
  expect_equal(unname(as.matrix(logfc(phip_obj))), replacement_matrix)
  expect_equal(unname(as.matrix(prob(phip_obj))), replacement_matrix)

  # Check that invalid replacement generates warning
  expect_error(counts(phip_obj) <- matrix(-1L, nrow = n_peptides, ncol = n_samples),
               "'counts' cannot have negative entries.")

  # Check peptideInfo replacement; new info has different rownames
  peptideInfo(phip_obj) <- virscan_info[, 1:10]
  pep_tidied <- .tidyPeptideInfo(virscan_info[, 1:10],
                                  dimnames(phip_obj)[[1]])
  pep_expect <- GenomicRanges::GRanges(seqnames = dimnames(phip_obj)[[1]],
     ranges = IRanges::IRanges(start = pep_tidied[["pep_start"]],
                               end = pep_tidied[["pep_end"]]))

  mcols(pep_expect) <- pep_tidied[["pep_meta"]]

  expect_equal(unname(peptideInfo(phip_obj)), pep_expect)

  # Check sampleInfo replacement, new info has different rownames
  sampleInfo$ART <- sample(c("yes", "no"), n_samples, replace = TRUE)
  sampleInfo(phip_obj) <- sampleInfo

  expect_equal(sampleInfo(phip_obj), sampleInfo)

})

test_that("assays can be added and removed", {

  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo, peptideInfo = virscan_info)

  replacement_matrix <- matrix(runif(n_samples*n_peptides, min = 1, max = 1e6),
                               nrow = n_peptides)

  expect_error(assay(phip_obj, "counts") <- NULL,
               paste0("`counts`, `logfc`, and `prob` assays must be included ",
                      "in a PhIPData object. The following assays are ",
                      "missing: counts"))
  assay(phip_obj, "new_assay") <- replacement_matrix
  expect_equal(unname(as.matrix(assay(phip_obj, "new_assay"))),
               replacement_matrix)
  assay(phip_obj, "new_assay") <- NULL
  expect_true(!"new_assay" %in% assayNames(phip_obj))

})

# Test coercion functions -----------
test_that("coercion functions work as expected", {
  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo, peptideInfo = virscan_info)

  expect_is(as(phip_obj, "list"), "list")
  expect_is(as(phip_obj, "List"), "List")
  expect_is(as(phip_obj, "DGEList"), "DGEList")

  expect_is(as(as(phip_obj, "list"), "PhIPData"), "PhIPData")
  expect_is(as(as(phip_obj, "List"), "PhIPData"), "PhIPData")
  expect_is(as(as(phip_obj, "DGEList"), "PhIPData"), "PhIPData")

})
