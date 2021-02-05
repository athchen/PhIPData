
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

The `PhIPData` class is used to store experimental results from
phage-immunoprecipitation sequencing (PhIP-set) experiments in a
matrix-like container.

Building on the
[`RangedSummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
class, `PhIPData` contains all of functionality of
`SummarizedExperiments` and includes additional operations to facilitate
analysis with PhIP-seq data. Like `SummarizedExperiments`, a key feature
of `PhIPData` is the coordination of metadata when subsetting `PhIPData`
objects. For example, if you wanted to examine experimental data for
peptides from one particular virus, you can subset the experimental data
and the associated peptide annotation with one command. This ensures all
metadata (for samples, peptides, etc.) remain synced with the
experimental data throughout analysis.

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
library(PhIPData)
```

## Components of a PhIPData Object

As reflected in the figure below, the structure of a `PhIPData` object
is nearly identical to the structure of a
`SummarizedExperiment`/`RangedSummarizedExperiment` object.

Each object contains at least three assays of data. These assays are:

  - `counts`: matrix of raw read counts,
  - `logfc`: matrix of log10 estimated fold-changes (in comparison to
    negative control samples),
  - `prob`: matrix of probabilities (p-values or posterior
    probabilities) associated with whether a sample shows an enriched
    antibody response to the particular peptide.

The rows of a `PhIPData` object represent peptides of interest and the
columns represent samples. Sample and peptide metadata are stored in
`DataFrame`s. Each row of the metadata `DataFrame` specifies the
peptide/sample, and the columns represent different features associated
with the peptides/samples.

In addition to sample- and peptide-specific metadata, experimental
metadata such as associated papers, experimental parameters, sequencing
dates, etc. are stored in a list-like component named `metadata`.

<div class="figure">

<img src="vignettes/extras/PhIPData.png" alt="Schematic of a PhIPData object. Commands used to access each component of the object are listed underneath its visual representation. Code in black indicates functions specific to `PhIPData` objects while functions in red extend `SummarizedExperiment` functions. Here, `pd` is a generic `PhIPData` object." width="\maxwidth" />

<p class="caption">

Schematic of a PhIPData object. Commands used to access each component
of the object are listed underneath its visual representation. Code in
black indicates functions specific to `PhIPData` objects while functions
in red extend `SummarizedExperiment` functions. Here, `pd` is a generic
`PhIPData` object.

</p>

</div>

## Example Use

For an example of how to use `PhIPData` objects, please see the package
vignette.
