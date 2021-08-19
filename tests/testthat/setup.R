## Simulated data to test PhIPData functions
set.seed(20210106)

## Import library
virscan_file <- system.file("extdata", "virscan.tsv", package = "PhIPData")
virscan_info <- readr::read_tsv(virscan_file,
    col_types = readr::cols(
        pep_id = readr::col_character(),
        pep_dna = readr::col_character(),
        pep_aa = readr::col_character(),
        pep_pos = readr::col_character(),
        pro_len = readr::col_double(),
        uniprot_acc = readr::col_character(),
        refseq = readr::col_character(),
        species = readr::col_character(),
        interspecies_specific = readr::col_character(),
        product = readr::col_character(),
        description = readr::col_character(),
        go = readr::col_character(),
        kegg = readr::col_character(),
        pfam = readr::col_character(),
        embl = readr::col_character(),
        interpro = readr::col_character(),
        pep_name = readr::col_character()
    )
) %>%
    as.data.frame()

n_samples <- 96L
n_peptides <- nrow(virscan_info)

## Generate data
counts <- matrix(sample(1:1e6, n_samples * n_peptides, replace = TRUE),
    nrow = n_peptides
)
logfc <- matrix(rnorm(n_samples * n_peptides, mean = 0, sd = 10),
    nrow = n_peptides
)
prob <- matrix(rbeta(n_samples * n_peptides, shape1 = 1, shape2 = 1),
    nrow = n_peptides
)

sampleInfo <- S4Vectors::DataFrame(
    sample_name = paste0("sample", 1:n_samples),
    gender = sample(c("M", "F"), n_samples,
        replace = TRUE
    ),
    group = sample(c("beads", "ART", "no ART"),
        n_samples,
        replace = TRUE
    )
)
## Set names
rownames(counts) <- rownames(logfc) <-
    rownames(prob) <- rownames(virscan_info) <-
    paste0("pep_", 1:n_peptides)

colnames(counts) <- colnames(logfc) <-
    colnames(prob) <- rownames(sampleInfo) <-
    paste0("sample_", 1:n_samples)

## default object
phip_obj <- PhIPData(
    counts = counts, logfc = logfc, prob = prob,
    sampleInfo = sampleInfo, peptideInfo = virscan_info
)

withr::defer(
    unlink(system.file(package = "PhIPData", "extdata/defaults.rda")),
    teardown_env()
)
