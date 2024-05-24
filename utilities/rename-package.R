#!/usr/bin/env -S Rscript --vanilla

# This script takes in a name and renames it to include the current R
# version and the current platform. It is intended to be used to rename
# packages that are being uploaded to releases on GitHub.

##  Parse the first argument as the package name + package version
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
    stop("Usage: rename-package.R <package-name>")
}

message("Received package name: ", args[[1L]])

packageName <- strsplit(args, "_")[[1L]][[1L]]
message("Package name: ", packageName)
packageVersion <- strsplit(args, "_")[[1L]][[2L]]
message("Package version: ", packageVersion)

newName <- paste(
    packageName, paste("R",
        getRversion(),
        sep = "-"
    ),
    paste0("v", packageVersion),
    sep = "_"
)

message("Renaming package to: ", newName)

value <- file.rename(args[[1L]], newName)

if (value) {
    message("Successfully renamed package to: ", newName)
    cat(newName)
} else {
    stop("Failed to rename package to: ", newName)
    cat("")
}
