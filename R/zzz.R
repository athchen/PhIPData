pkg_env <- new.env()

.createDefaults <- function() {

    # Create default variable for beads name
    beads_path <- BiocFileCache::bfcnew(
        pkg_env$beer_cache, "beads_name",
        ext = ".rds"
    )
    saveRDS("beads", beads_path)

    # Copy default alias list to cache
    BiocFileCache::bfcadd(pkg_env$beer_cache, "alias",
        system.file(package = "PhIPData", "extdata/alias.rda"),
        action = "copy"
    )

    # Copy small VirScan library to cache
    BiocFileCache::bfcadd(pkg_env$beer_cache, "library_virscan",
        system.file(package = "PhIPData", "libraries", "virscan.rds"),
        action = "copy"
    )
    TRUE
}

.onLoad <- function(libname, pkgname) {

    ## Load cache
    cache_path <- tools::R_user_dir("beer", which = "cache")
    assign("beer_cache",
        BiocFileCache::BiocFileCache(cache_path, ask = FALSE),
        envir = pkg_env
    )

    ## If cache does not exist, add defaults
    if (length(pkg_env$beer_cache) == 0) {
        tmp <- .createDefaults()
    }

    ## Load beads name
    beads_path <- BiocFileCache::bfcquery(pkg_env$beer_cache, "beads_name")$rpath
    pkg_env$BEADS_NAME <- readRDS(beads_path)

    ## Load alias
    alias_path <- BiocFileCache::bfcquery(pkg_env$beer_cache, "alias")$rpath
    load(alias_path, envir = pkg_env)
}
