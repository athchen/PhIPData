context("Libraries can be stored and loaded for populating peptide information.")

test_that("libraries can be created and used to make valid PhIPData objects.", {
  is_windows <- grepl("windows", .Platform$OS.type)

  library_loc <- system.file("libraries", package = "PhIPData")
  library_loc <- if(is_windows){
    gsub("/", "\\", library_loc, fixed = TRUE)
    } else { library_loc }
  expect_equal(getLibraryPath(), library_loc)

  # Test library path functions
  expect_error(setLibraryPath("invalid_path"), "Invalid specified path.")

  extdata_loc <- system.file("extdata", package = "PhIPData")
  setLibraryPath(extdata_loc)
  expect_equal(getLibraryPath(), extdata_loc)

  # Test library creation
  if(file.exists(paste0(extdata_loc, "/virscan.rds"))){
    file.remove(paste0(extdata_loc, "/virscan.rds"))
  }

  virscan_file <- system.file("extdata", "virscan.tsv", package = "PhIPData")
  virscan_info <- readr::read_tsv(virscan_file,
                                  col_types = readr::cols(
                                    pep_id = readr::col_character(),
                                    pro_id = readr::col_character(),
                                    pos_start = readr::col_double(),
                                    pos_end = readr::col_double(),
                                    UniProt_acc = readr::col_character(),
                                    pep_dna = readr::col_character(),
                                    pep_aa = readr::col_character(),
                                    pro_len = readr::col_double(),
                                    taxon_id = readr::col_double(),
                                    species = readr::col_character(),
                                    genus = readr::col_character(),
                                    product = readr::col_character()
                                  )) %>%
    as.data.frame()

  makeLibrary(virscan_info, "virscan")
  expect_true(file.exists(paste0(extdata_loc, "/virscan.rds")))

  # test use function
  n_samples <- 96L
  n_peptides <- nrow(virscan_info)
  counts <- matrix(sample(1:1e6, n_samples*n_peptides, replace = TRUE),
                   nrow = n_peptides)
  logfc <- matrix(rnorm(n_samples*n_peptides, mean = 0, sd = 10),
                  nrow = n_peptides)
  prob <- matrix(rbeta(n_samples*n_peptides, shape1 = 1, shape2 = 1),
                 nrow = n_peptides)

  sampleInfo <- DataFrame(sample_name = paste0("sample", 1:n_samples),
                          gender = sample(c("M", "F"), n_samples,
                                          replace = TRUE))

  rownames(counts) <- rownames(logfc) <-
    rownames(prob) <- rownames(virscan_info) <-
    paste0("pep_", 1:n_peptides)

  colnames(counts) <- colnames(logfc) <-
    colnames(prob) <- rownames(sampleInfo) <-
    paste0("sample_", 1:n_samples)

  phip_obj <- PhIPData(counts = counts, logfc = logfc, prob = prob,
                       sampleInfo = sampleInfo,
                       peptideInfo = getLibrary("virscan"))

  # clean-up test space
  Sys.unsetenv("PHIP_LIBRARY_PATH")
  file.remove(paste0(extdata_loc, "/virscan.rds"))
})
