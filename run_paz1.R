#!/usr/bin/env Rscript

# Load the pipeline function
source("R/koios_pipeline.R")

cat("=== Processing Paz1 VCF with Preprocessing ===\n")
cat("  - Input directory: input/\n")
cat("  - Output directory: output/\n")
cat("  - hg19 -> hg38 liftover\n")
cat("  - Quality filtering (AF >= 1%, DP >= 100)\n")
cat("  - Remove undefined genotypes (./.)  \n")
cat("  - Separate CNV from SNP/INDEL\n")
cat("  - Select most representative allele\n\n")

# Run pipeline on Paz1 VCF
koios_vrs_pipeline(
    vcf_path = "Paz1_v1_Non-Filtered_2025-10-30_04_31_11.vcf",  # Will look in input/ directory
    output_base_name = "Paz1_output",
    ref_genome = "hg19",  # Input VCF is hg19, will perform automatic liftover
    af_threshold = 0.01,  # 1% allele frequency threshold for gnomAD filtering
    min_dp = 100,  # Minimum read depth for quality filtering
    input_dir = "input",  # Input directory containing VCF and CSV
    output_dir = "output",  # Output directory for all generated files
    default_vus_id = 1028197L  # OMOP concept ID for VUS
)

cat("\n=== Pipeline completed successfully ===\n")
cat("Output files generated:\n")
cat("  - Paz1_output_VUS.csv (Variants of Unknown Significance)\n")
cat("  - Paz1_output_Note.csv (Clinically significant variants)\n")
cat("  - Paz1_output_Phenopackets.json (GA4GH Phenopackets v2)\n")
cat("  - Paz1_v1_Non-Filtered_2025-10-30_04_31_11_hg38.vcf (Lifted VCF)\n")
