# PhIPData

<!-- badges: start -->
[![R build status](https://github.com/athchen/PhIPData/workflows/R-CMD-check/badge.svg)](https://github.com/athchen/PhIPData/actions)
<!-- badges: end -->

# To Do
* Make a show method
  + What do we want to include in the show method? Is the default okay?
* Make subset functions
  + Make function to subset peptides by virus
* Make library functions
* Add documentation
  + Note that many `SummarizedExperiment` getters and setters will also work on the PhIPData objects
      - `metadata()`: returns metadata
      - `dimnames()`: returns sample names and peptide names
      - `rownames()`: sets all the sample names
      - `colnames()`: sets all the peptide names
* Make functions to select default libraries (e.g they can pass `use_library("virscan")`)
* Make function to get experiment summary (e.g. library sizes)
* `pkgdown` template edits
