#!/usr/bin/env Rscript
#
# create_clinical_csv.R
#
# KOIOS-VRS Pipeline: Clinical Metadata CSV Helper
#
# Authors:
#   - Gabriele Bucci (Lead Developer, Project Ideator)
#   - Michela Riba (Project Ideator)
#   - Guido Scicolone (Developer, Project Ideator)
#   - Muhammad Tajwar (Developer)
#
# Part of: KOIOS-VRS Pipeline (https://github.com/gbucci/koios_vrs)
#
# Description:
#   Helper script to create clinical metadata CSV files for the KOIOS-VRS pipeline.
#   Can be run in two modes:
#   1. Interactive: Run without arguments, will prompt for each field
#   2. Command-line: Pass arguments directly
#
################################################################################

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    # Interactive mode
    cat("\n=== CREATE CLINICAL METADATA CSV ===\n\n")

    # Get VCF path
    vcf_path <- readline(prompt = "Enter VCF file path (e.g., 'input/sample.vcf'): ")
    if (!file.exists(vcf_path)) {
        stop(sprintf("VCF file not found: %s", vcf_path))
    }

    # Determine CSV path
    csv_path <- sub("\\.vcf$", ".csv", vcf_path, ignore.case = TRUE)
    cat(sprintf("CSV will be created at: %s\n\n", csv_path))

    # Get sample ID
    vcf_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
    sample_id <- colnames(vcf_obj@gt)[2]
    cat(sprintf("Detected sample ID from VCF: %s\n\n", sample_id))

    # Get clinical data
    diagnosis_label <- readline(prompt = "Diagnosis label (e.g., 'Lung Cancer'): ")
    diagnosis_ncit_id <- readline(prompt = "Diagnosis NCIT ID (e.g., 'NCIT:C7377'): ")

    cat("\nSex options:\n")
    cat("  1. MALE (PATO:0000384)\n")
    cat("  2. FEMALE (PATO:0000383)\n")
    sex_choice <- readline(prompt = "Select (1 or 2): ")

    if (sex_choice == "1") {
        sex_label <- "MALE"
        sex_pato_id <- "PATO:0000384"
    } else if (sex_choice == "2") {
        sex_label <- "FEMALE"
        sex_pato_id <- "PATO:0000383"
    } else {
        stop("Invalid sex choice")
    }

    age_of_onset <- readline(prompt = "\nAge of onset in years (e.g., '45'): ")
    age_of_onset_iso <- sprintf("P%sY", age_of_onset)

} else {
    # Command-line mode
    if (length(args) < 6) {
        cat("Usage: Rscript create_clinical_csv.R <vcf_path> <diagnosis_label> <diagnosis_ncit_id> <sex_label> <sex_pato_id> <age_years>\n")
        cat("\nExample:\n")
        cat("  Rscript create_clinical_csv.R input/sample.vcf 'Lung Cancer' NCIT:C7377 MALE PATO:0000384 45\n")
        stop("Insufficient arguments")
    }

    vcf_path <- args[1]
    if (!file.exists(vcf_path)) {
        stop(sprintf("VCF file not found: %s", vcf_path))
    }

    csv_path <- sub("\\.vcf$", ".csv", vcf_path, ignore.case = TRUE)

    # Get sample ID from VCF
    vcf_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
    sample_id <- colnames(vcf_obj@gt)[2]

    diagnosis_label <- args[2]
    diagnosis_ncit_id <- args[3]
    sex_label <- args[4]
    sex_pato_id <- args[5]
    age_of_onset_iso <- sprintf("P%sY", args[6])
}

# Create data frame
clinical_data <- data.frame(
    sample_id = sample_id,
    diagnosis_label = diagnosis_label,
    diagnosis_ncit_id = diagnosis_ncit_id,
    sex_label = sex_label,
    sex_pato_id = sex_pato_id,
    age_of_onset_iso = age_of_onset_iso,
    stringsAsFactors = FALSE
)

# Show summary
cat("\n=== CLINICAL METADATA SUMMARY ===\n")
print(clinical_data)
cat("\n")

# Write CSV
write.csv(clinical_data, file = csv_path, row.names = FALSE, quote = TRUE)
cat(sprintf("✓ CSV file created: %s\n\n", csv_path))

cat("You can now run the pipeline:\n")
cat(sprintf("  Rscript run_pipeline.R %s output_name\n\n", vcf_path))
