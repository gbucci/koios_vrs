#!/usr/bin/env Rscript

# Load the pipeline function
source("R/koios_pipeline.R")

cat("=== Processing Melanoma Example ===\n")
cat("  - Input directory: examples/\n")
cat("  - Output directory: output/\n")
cat("  - Reference genome: hg38 (no liftover needed)\n")
cat("  - Quality filtering (AF >= 1%, DP >= 50) - relaxed for example\n")
cat("  - Remove undefined genotypes (./.)  \n")
cat("  - Separate CNV from SNP/INDEL\n")
cat("  - Select most representative allele\n\n")

# Run pipeline on synthetic melanoma sample
# Note: Using relaxed thresholds (min_dp = 50) for this small example dataset
koios_vrs_pipeline(
    vcf_path = "melanoma_sample.vcf",
    output_base_name = "melanoma_output",
    ref_genome = "hg38",  # Input VCF is already hg38
    af_threshold = 0.01,  # 1% allele frequency threshold for gnomAD filtering
    min_dp = 50,  # Minimum read depth (relaxed to 50 for example data)
    input_dir = "examples",  # Input directory containing example VCF and CSV
    output_dir = "output",  # Output directory for all generated files
    default_vus_id = 1028197L  # OMOP concept ID for VUS
)

cat("\n=== Pipeline completed successfully ===\n")
cat("Output files generated in output/ directory:\n")
cat("  - melanoma_output_VUS.csv (Variants of Unknown Significance)\n")
cat("  - melanoma_output_Note.csv (Clinically significant variants)\n")
cat("  - melanoma_output_Phenopackets.json (GA4GH Phenopackets v2)\n")
cat("  - melanoma_output_GeneSummary.csv (Gene-level mutation summary)\n")
cat("  - melanoma_output_annotated.vcf (Final annotated VCF)\n")
cat("  - melanoma_sample_preprocessed.vcf (Quality-filtered VCF)\n")
cat("  - melanoma_sample_CNV.vcf (Separated CNV variants)\n")
