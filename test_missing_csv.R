#!/usr/bin/env Rscript

# Test: What happens when CSV is missing in Rscript mode

library(vcfR)
library(KOIOS)

source("R/koios_pipeline.R")

cat("\n=== TESTING MISSING CSV BEHAVIOR ===\n\n")

# Temporarily rename the CSV to simulate missing file
csv_file <- "RWD/Paz1_v1_Non-Filtered_2025-10-30_04_31_11.csv"
csv_backup <- "RWD/Paz1_v1_Non-Filtered_2025-10-30_04_31_11.csv.backup"

if (file.exists(csv_file)) {
    file.rename(csv_file, csv_backup)
    cat("Temporarily hidden CSV file\n\n")
}

# Try to run pipeline - should fail with helpful message
tryCatch({
    result <- koios_vrs_pipeline(
        vcf_path = "RWD/Paz1_v1_Non-Filtered_2025-10-30_04_31_11.vcf",
        output_base_name = "test_missing_csv",
        ref_genome = "hg19"
    )
}, error = function(e) {
    cat("\n=== EXPECTED ERROR ===\n")
    cat("Error message:\n")
    cat(conditionMessage(e), "\n")
})

# Restore CSV
if (file.exists(csv_backup)) {
    file.rename(csv_backup, csv_file)
    cat("\n\nRestored CSV file\n")
}

cat("\n=== TEST COMPLETE ===\n")
