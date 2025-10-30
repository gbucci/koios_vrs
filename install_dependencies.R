#!/usr/bin/env Rscript

cat("Installing required Bioconductor packages for hg19 liftover support...\n\n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install required Bioconductor packages
packages <- c("rtracklayer", "GenomicRanges", "AnnotationHub", "IRanges")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("Installing %s...\n", pkg))
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
        cat(sprintf("%s is already installed.\n", pkg))
    }
}

cat("\n=== All dependencies installed successfully ===\n")
