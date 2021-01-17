context("Create, modify, and delete aliases for peptide subsetting")

test_that("the alias database can be accessed and modified", {

  skip("Reworking class")

  # Test path functions
  expect_error(setAliasPath("bad_path"), "Invalid alias file location.")
  expect_error(setAliasPath(system.file("extdata/virscan.tsv",
                                        package = "PhIPData")),
               "Invalid file type.")

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
  expect_error(getAlias("test_virus"),
               "Virus does not exist in the alias database.")
  expect_error(deleteAlias("test_virus"),
               "Virus does not exist in the alias database.")

})
