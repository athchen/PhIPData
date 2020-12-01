context("PhIPData constructor works for a variety of poorly-defined inputs")

test_that("a PhIPData object can be created if one or more of the parameters are missing.", {

  # Set up information for test cases
  # virscan_info <- readr::read_tsv("tests/testdata/VirScan_annotation.tsv",
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

  rownames(counts) <- rownames(logfc) <- rownames(prob) <- paste0("pep_", 1:n_peptides)
  colnames(counts) <- colnames(logfc) <- colnames(prob) <- paste0("sample_", 1:n_samples)

  sampleInfo <- DataFrame(sample_name = paste0("sample", 1:n_samples),
                          gender = sample(c("M", "F"), n_samples, replace = TRUE))

  expect_is(PhIPData(), "PhIPData")
  expect_is(PhIPData(counts = counts), "PhIPData")
  expect_is(PhIPData(logfc = logfc), "PhIPData")
  expect_is(PhIPData(prob = prob), "PhIPData")
  expect_is(PhIPData(counts = counts, logfc = logfc), "PhIPData")
  expect_is(PhIPData(counts = counts, prob = prob), "PhIPData")
  expect_is(PhIPData(logfc = logfc, prob = prob), "PhIPData")
  expect_error(PhIPData(peptideInfo = virscan_info), "Cannot create empty PhIPData object with only one of sampleInfo or peptideInfo.")
  expect_error(PhIPData(sampleInfo = sampleInfo), "Cannot create empty PhIPData object with only one of sampleInfo or peptideInfo.")
})


