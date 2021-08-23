pkg_env <- new.env()

.createDefaults <- function() {
    assign("BEADS_NAME", "beads", envir = pkg_env)
    assign("ALIAS_PATH",
        system.file(package = "PhIPData", "extdata/alias.rda"),
        envir = pkg_env
    )
    assign("PHIP_LIBRARY_PATH",
        system.file(package = "PhIPData", "libraries"),
        envir = pkg_env
    )
    save(
        list = c("BEADS_NAME", "ALIAS_PATH", "PHIP_LIBRARY_PATH"),
        file = paste0(
            system.file(package = "PhIPData", "extdata"),
            "/defaults.rda"
        ),
        envir = pkg_env
    )
}

.onLoad <- function(libname, pkgname) {
    # Build default variables by system if file does not exists
    extdata_path <- system.file(package = "PhIPData", "extdata")
    if (system.file(package = "PhIPData", "extdata/defaults.rda") == "") {
        .createDefaults()
    }
    load(system.file(package = "PhIPData", "extdata/defaults.rda"),
        envir = pkg_env
    )

    ## Check if defaults info is valid, if not, rebuild with defaults
    if (any(!file.exists(c(pkg_env$ALIAS_PATH, pkg_env$PHIP_LIBRARY_PATH)))) {
        .createDefaults()
    }
    load(get("ALIAS_PATH", envir = pkg_env), envir = pkg_env)
}
