
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PhIPData

<!-- badges: start -->

[![R build
status](https://github.com/athchen/PhIPData/workflows/R-CMD-check/badge.svg)](https://github.com/athchen/PhIPData/actions)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/PhIPData.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/PhIPData)
[![Codecov test
coverage](https://codecov.io/gh/athchen/PhIPData/branch/main/graph/badge.svg)](https://codecov.io/gh/athchen/PhIPData?branch=main)
<!-- badges: end -->

The goal of PhIPData is to …

# To Do

  - Make a show method
      - What do we want to include in the show method? Is the default
        okay?
  - Make subset functions
      - Make function to subset peptides by virus
  - Make library functions
  - Add documentation
      - Note that many `SummarizedExperiment` getters and setters will
        also work on the PhIPData objects
          - `metadata()`: returns metadata
          - `dimnames()`: returns sample names and peptide names
          - `rownames()`: sets all the sample names
          - `colnames()`: sets all the peptide names
  - Make functions to select default libraries (e.g they can pass
    `use_library("virscan")`)
  - Make function to get experiment summary (e.g. library sizes)
  - `pkgdown` template edits

## Installation

We recommend installing the stable release version of `PhIPData` in
Bioconductor. This can be done using `BiocManager`:

``` r
if (!require("BiocManager"))
    install.packages("BiocManager")
    
BiocManager::install("PhIPData")
```

To load the package:

``` r
#library(PhIPData)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(PhIPData)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
