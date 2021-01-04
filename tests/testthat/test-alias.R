context("Create, modify, and delete aliases for peptide subsetting")

test_that("the alias database can be accessed and modified", {

  # check setter when virus does not exist
  setAlias("test_virus", "test_pattern")
  expect_equal(getAlias("test_virus"), "test_pattern")

  # Check setter when virus exists
  setAlias("test_virus", "test_pattern2")
  expect_equal(getAlias("test_virus"), "test_pattern2")

  # Check remove function
  deleteAlias("test_virus")
  expect_error(getAlias("test_virus"), "Virus does not exist in alias database.")

})
