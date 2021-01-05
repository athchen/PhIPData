context("Library creators and accessors work as expected.")

test_that("libraries can be created and stored to make valid PhIPData objects.", {
  expect_error(getLibraryPath(), "Library path does not exist.")

  setLibraryPath("../testdata/")
  expect_equal(getLibraryPath(), normalizePath("../testdata/"))

  virscan_info <- readr::read_tsv("../testdata/VirScan_annotation.tsv",
                                  col_types = readr::cols(
                                    .default = readr::col_character(),
                                    pep_id = readr::col_double(),
                                    pep_rank = readr::col_double(),
                                    pos_start = readr::col_double(),
                                    pos_end = readr::col_double(),
                                    pro_len = readr::col_double()
                                  )) %>%
    as.data.frame()

  makeLibrary(virscan_info, "virscan")
  expect_true(file.exists("../testdata/virscan.rds"))

  # test use function
  n_samples <- 96L
  n_peptides <- nrow(virscan_info)
  counts <- matrix(runif(n_samples*n_peptides, min = 1, max = 1e6),
                   nrow = n_peptides)
  logfc <- matrix(rnorm(n_samples*n_peptides, mean = 0, sd = 10),
                  nrow = n_peptides)
  prob <- matrix(rbeta(n_samples*n_peptides, shape1 = 1, shape2 = 1),
                 nrow = n_peptides)

  sampleInfo <- DataFrame(sample_name = paste0("sample", 1:n_samples),
                          gender = sample(c("M", "F"), n_samples, replace = TRUE))

  rownames(counts) <- rownames(logfc) <- rownames(prob) <- rownames(virscan_info) <- paste0("pep_", 1:n_peptides)
  colnames(counts) <- colnames(logfc) <- colnames(prob) <- rownames(sampleInfo) <- paste0("sample_", 1:n_samples)

  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo, peptideInfo = getLibrary("virscan"))
  # clean-up test space
  Sys.setenv(LIBRARY_PATH = "")
})
