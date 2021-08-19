context("Check PhIPData subsetting wrappers")

# source("setup.R")

test_that("beads-only subsetting works", {
    # empty case
    expect_is(subsetBeads(PhIPData()), "PhIPData")

    # full case
    phip_obj <- PhIPData(
        counts = counts, logfc = logfc, prob = prob,
        sampleInfo = sampleInfo, peptideInfo = virscan_info
    )
    beads_sub <- subsetBeads(phip_obj)
    expect_equal(
        dim(beads_sub),
        c(nrow(phip_obj), sum(sampleInfo$group == "beads"))
    )
    expect_equal(rownames(phip_obj), rownames(beads_sub))
    expect_equal(
        colnames(phip_obj)[sampleInfo$group == "beads"],
        colnames(beads_sub)
    )
})
