context("PhIPData constructor works for a variety of poorly-defined inputs")

virscan_info <- readr::read_tsv("../testdata/VirScan_annotation.tsv",
                                col_types = readr::cols(
                                  .default = readr::col_character(),
                                  pep_id = readr::col_double(),
                                  pep_rank = readr::col_double(),
                                  pos_start = readr::col_double(),
                                  pos_end = readr::col_double(),
                                  pro_len = readr::col_double()
                                ))
n_samples <- 96L
n_peptides <- nrow(virscan_info)
counts <- matrix(runif(n_samples*n_peptides, min = 1, max = 1e6), nrow = n_peptides)
logfc <- matrix(rnorm(n_samples*n_peptides, mean = 0, sd = 10), nrow = n_peptides)
prob <- matrix(rbeta(n_samples*n_peptides, shape1 = 1, shape2 = 1), nrow = n_peptides)

sampleInfo <- DataFrame(sample_name = paste0("sample", 1:n_samples),
                        gender = sample(c("M", "F"), n_samples, replace = TRUE))

rownames(counts) <- rownames(logfc) <- rownames(prob) <- rownames(virscan_info) <- paste0("pep_", 1:n_peptides)
colnames(counts) <- colnames(logfc) <- colnames(prob) <- rownames(sampleInfo) <- paste0("sample_", 1:n_samples)

test_that("a PhIPData object can be created if one or more of the parameters are missing.", {

  # Set up information for test cases
  # virscan_info <- readr::read_tsv("tests/testdata/VirScan_annotation.tsv",

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
  expect_is(PhIPData(counts = counts, peptideInfo = virscan_info, sampleInfo = sampleInfo), "PhIPData")

  # only two assays are present
  expect_is(PhIPData(counts = counts, logfc = logfc), "PhIPData")
  expect_is(PhIPData(counts = counts, prob = prob), "PhIPData")
  expect_is(PhIPData(logfc = logfc, prob = prob), "PhIPData")

  # only two assay and metadata
  expect_is(PhIPData(counts = counts, logfc = logfc, sampleInfo = sampleInfo), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc, peptideInfo = virscan_info), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc,
                     peptideInfo = virscan_info, sampleInfo = sampleInfo), "PhIPData")

  # all 3 assays
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob), "PhIPData")

  # all 3 assays and metadata
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     sampleInfo = sampleInfo), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     peptideInfo = virscan_info), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc, prob = prob,
                     peptideInfo = virscan_info, sampleInfo = sampleInfo), "PhIPData")
})

test_that("getter functions operate as expected.", {
  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo, peptideInfo = virscan_info)

  expect_is(counts(phip_obj), "DataFrame")
  expect_is(logfc(phip_obj), "DataFrame")
  expect_is(prob(phip_obj), "DataFrame")
  expect_is(peptideInfo(phip_obj), "GRanges")
  expect_is(sampleInfo(phip_obj), "DataFrame")
})


