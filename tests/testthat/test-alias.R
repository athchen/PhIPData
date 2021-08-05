context("Aliases for peptide subsetting can be created, modified, and deleted.")

test_that("the alias database can be accessed and modified", {

    # Test path functions
    expect_equal(
        getAliasPath(),
        system.file(package = "PhIPData", "extdata/alias.rda")
    )
    expect_error(setAliasPath("bad_path"), "Invalid alias file location.")
    expect_error(
        setAliasPath(system.file("extdata/virscan.tsv",
            package = "PhIPData"
        )),
        "Invalid file type."
    )

    alias_path <- system.file(package = "PhIPData", "extdata/alias.rda")
    setAliasPath(alias_path)
    expect_equal(Sys.getenv("ALIAS_PATH"), alias_path)
    expect_equal(getAliasPath(), alias_path)

    # check setter when virus does not exist
    setAlias("test_virus", "test_pattern")
    expect_equal(getAlias("test_virus"), "test_pattern")

    # Check setter when virus exists
    setAlias("test_virus", "test_pattern2")
    expect_equal(getAlias("test_virus"), "test_pattern2")

    # Check remove function
    deleteAlias("test_virus")
    expect_equal(getAlias("test_virus"), NA_character_)
    expect_error(
        deleteAlias("test_virus"),
        "Virus does not exist in the alias database."
    )

    # clean-up test space
    Sys.unsetenv("ALIAS_PATH")
})

test_that("Alias functions can be applied to vectors", {
    # Check vectorization
    expect_equal(getAlias(c("HIV", "EBV")), c("Human immunodeficiency virus", "Epstein-Barr"))
    setAlias(c("test_1", "HIV"), c("pattern_1", "hiv"))
    expect_equal(getAlias(c("test_1", "HIV")), c("pattern_1", "hiv"))
    deleteAlias(c("test_1", "HIV"))
    expect_equal(getAlias(c("test_1", "HIV")), c(NA_character_, NA_character_))

    ## clean-up test space
    setAlias("HIV", "Human immunodeficiency virus")
})
