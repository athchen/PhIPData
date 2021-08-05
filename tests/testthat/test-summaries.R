context("Computing summaries of PhIPData objects")

source("setup.R")

test_that("library sizes are accurately calculated", {
    expect_equal(librarySize(phip_obj), colSums(counts))
    expect_equal(
        librarySize(phip_obj, withDimnames = FALSE),
        unname(colSums(counts))
    )
})

test_that("proportion of reads pulled are accurately calculated", {
    manual_prop <- counts / matrix(rep(colSums(counts), nrow(counts)),
        nrow = nrow(counts), byrow = TRUE
    )
    expect_equal(propReads(phip_obj), manual_prop)
    expect_equal(
        propReads(phip_obj, withDimnames = FALSE),
        unname(manual_prop)
    )
})
