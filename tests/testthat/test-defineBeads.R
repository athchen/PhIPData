context("Varying strings can be used to indicate beads-only samples")

test_that("varying strings can be used to indicate beads-only samples", {
    expect_equal(getBeadsName(), "beads")

    # Test setters
    setBeadsName("BEADS_ONLY")
    expect_equal(getBeadsName(), "BEADS_ONLY")
    setBeadsName(1)
    expect_equal(getBeadsName(), "1")
    expect_error(setBeadsName(NA), "Beads cannot be specified via NA.")

    # Clean environment
    setBeadsName("beads")
})
