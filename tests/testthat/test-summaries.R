context("Computing summaries of PhIPData objects")

source("setup.R")

test_that("library sizes are accurately calculated", {
  expect_equal(librarySize(phip_obj), colSums(counts))
  expect_equal(librarySize(phip_obj, withDimnames = FALSE),
               unname(colSums(counts)))
})
