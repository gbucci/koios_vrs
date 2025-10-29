# run_example.R

# NOTE: This script ASSUMES it is being run from the project root directory
# (i.e., from the 'koiosVrsCancerPipe' folder).
# Execution command should be: Rscript inst/scripts/run_example.R

# Source the pipeline function from the R directory
source("R/koios_pipeline.R") 

# Define paths for the melanoma example relative to the project root
vcf_file <- "inst/extdata/melanoma_sample.vcf"
output_name <- "melanoma_test_output"

if (!file.exists(vcf_file)) {
    # If the VCF doesn't exist, it means the script was run from the wrong directory.
    stop(paste("Example VCF not found at expected path:", vcf_file, 
               "Please ensure you run this script from the root directory: koiosVrsCancerPipe/"))
}

# --- Execution ---

cat("--- Starting pipeline for:", vcf_file, "---\n")

# Run the pipeline function
results <- koios_vrs_pipeline(
    vcf_path = vcf_file,
    output_base_name = output_name,
    af_threshold = 0.01 
)

cat("\n--- Pipeline Complete ---\n")
cat("Output Files:\n")
cat("- Phenopackets JSON:", paste0(output_name, "_Phenopackets.json"), "\n")
cat("- VUS Variants CSV:", paste0(output_name, "_VUS.csv"), "\n")
cat("- Note Variants CSV:", paste0(output_name, "_Note.csv"), "\n")
