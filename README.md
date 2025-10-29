# koiosVrsCancerPipe

[![DOI](https://zenodo.org/badge/1072319376.svg)](https://doi.org/10.5281/zenodo.17466333)

## Overview

This R package implements a robust pipeline to process somatic or germline cancer variants from a VCF (Variant Call Format) file, enrich them with **KOIOS** (HGVSg, ClinGen, OMOP concepts), retrieve **VRS identifiers** from the external registry, filter out high-frequency benign variants, and generate interoperable output files (VUS/Note CSVs and **Phenopackets v2 JSON**).

The pipeline uses patient-specific clinical metadata provided in a companion CSV file to populate the Phenopackets.

### Key Acronyms

- **VCF**: Variant Call Format - standard file format for storing genomic variant calls
- **KOIOS**: Knowledge Organization for Interpretation of Omics Sequences - framework for variant enrichment
- **VRS**: Variant Representation Specification - standard for representing genomic variants
- **HGVSg**: Human Genome Variation Society genomic nomenclature
- **OMOP**: Observational Medical Outcomes Partnership - common data model
- **VUS**: Variant of Unknown Significance
- **AF**: Allele Frequency
- **NCIT**: National Cancer Institute Thesaurus
- **PATO**: Phenotype And Trait Ontology

## Package Structure

This is a standard R package following Bioconductor conventions:
- `R/koios_pipeline.R`: Main pipeline function (`koios_vrs_pipeline()`)
- `inst/extdata/`: Example input files (VCF and clinical CSV files for melanoma, breast, and CRC samples)
- `inst/scripts/run_example.R`: Example script demonstrating pipeline usage
- `man/`: Generated R documentation (roxygen2)
- `vignettes/`: Package documentation in Rmarkdown format
- `DESCRIPTION`: Package metadata and dependencies
- `NAMESPACE`: Exported functions (managed by roxygen2)

## Requirements

**R Environment:** R (>= 4.0), Bioconductor

**Required Packages:**
- `KOIOS`: Variant enrichment framework
- `httr`: HTTP requests to VRS API
- `jsonlite`: JSON parsing and generation
- `vcfR`: VCF file parsing
- `uuid`: UUID generation for Phenopackets
- `utils`, `stats`: Base R utilities

**Input Files:**
- VCF files must be in **hg38** reference genome format
- A companion clinical metadata CSV file must exist with the same base name as the VCF (e.g., `sample.vcf` requires `sample.csv`)
- VCF must contain exactly one sample column
- The sample ID in the VCF must match the `sample_id` in the clinical CSV

## Installation

```R
# Install required packages
install.packages(c("httr", "jsonlite", "vcfR", "uuid"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("KOIOS")

# Install from source
devtools::install()

# Or use standard R CMD INSTALL
R CMD INSTALL .
```

## Usage

### Run Example Test

This example demonstrates the pipeline using the provided sample data in `inst/extdata`:

```R
# Source the main pipeline function
source("R/koios_pipeline.R")

# Run example script
Rscript inst/scripts/run_example.R
```

### Using the Main Function

```R
source("R/koios_pipeline.R")

koios_vrs_pipeline(
    vcf_path = "path/to/sample.vcf",
    output_base_name = "output_prefix",
    af_threshold = 0.01,  # 1% allele frequency threshold (default)
    default_vus_id = 1028197L  # OMOP concept ID for VUS (default)
)
```

### Clinical Metadata CSV Format

Required columns:
- `sample_id`: Must match the sample column name in the VCF
- `diagnosis_label`: Disease name (e.g., "Melanoma")
- `diagnosis_ncit_id`: NCIT ontology ID (e.g., "NCIT:C3224")
- `sex_label`: "MALE" or "FEMALE"
- `sex_pato_id`: PATO ontology ID (e.g., "PATO:0000383")
- `age_of_onset_iso`: ISO 8601 duration format (e.g., "P60Y" for 60 years)

Example: `inst/extdata/melanoma_sample.csv`

## Pipeline Architecture

The `koios_vrs_pipeline()` function executes these sequential steps:

1. **Setup and Input Validation** (R/koios_pipeline.R:17-38)
   - Loads VCF and validates single-sample requirement
   - Locates and validates companion clinical CSV file
   - Extracts sample ID from VCF

2. **KOIOS Processing** (R/koios_pipeline.R:40-54)
   - Loads KOIOS concepts and hg38 reference genome
   - Processes VCF into dataframe format
   - Generates HGVSg nomenclature for variants
   - Enriches with ClinGen data
   - Adds OMOP concepts for variant classification
   - Assigns default VUS classification to unclassified variants

3. **Allele Frequency Filtering** (R/koios_pipeline.R:56-70)
   - Filters variants based on `koios_gnomAD_AF` population frequency
   - Default threshold: AF > 1% variants are excluded
   - High-frequency benign variants removed from downstream analysis

4. **VRS API Query** (R/koios_pipeline.R:72-124)
   - Queries VRS registry API (`https://reg.genome.network`) for each variant
   - Uses HGVSg notation to retrieve VRS allele IDs
   - Implements retry logic with exponential backoff for API failures
   - Falls back to alternative HGVSg column (`koios_hgvsg`) if primary lookup fails
   - Adds `vrs_id` and `vrs_seq_id` fields to variant dataframe

5. **Phenopackets v2 JSON Generation** (R/koios_pipeline.R:126-181)
   - Creates GA4GH Phenopackets v2 compliant JSON for each variant
   - Integrates clinical metadata (diagnosis, sex, age of onset)
   - Embeds VRS allele IDs and OMOP classifications
   - Includes genomic file references
   - Generates unique UUIDs for each phenopacket

6. **CSV Export** (R/koios_pipeline.R:183-189)
   - Separates variants into VUS and clinically significant (Note)
   - Exports two CSV files with complete variant annotations

## Output Files

For output base name `{prefix}`, the pipeline generates:
- `{prefix}_Phenopackets.json`: Array of Phenopackets v2 objects
- `{prefix}_VUS.csv`: Variants of Unknown Significance with all annotations
- `{prefix}_Note.csv`: Clinically significant variants (non-VUS)

## Development

### Generate Documentation

```R
# Generate documentation from roxygen2 comments
roxygen2::roxygenize()
```

### Testing with Different Samples

Three example datasets are provided in `inst/extdata/`:
- `melanoma_sample.{vcf,csv}`
- `breast_sample.{vcf,csv}`
- `crc_sample.{vcf,csv}` (colorectal cancer)

Modify `inst/scripts/run_example.R` to test different samples.

## API Dependencies

**VRS Registry API**: `https://reg.genome.network/vrAllele?hgvs={hgvs_string}`
- Rate limiting: 0.05s delay between requests (R/koios_pipeline.R:120)
- Retry strategy: Up to 3 attempts with exponential backoff for 429/5xx errors
- Timeout: 10 seconds per request

## Common Pitfalls

1. **Reference genome mismatch**: VCF must be hg38. hg19/GRCh37 VCFs will fail during KOIOS processing.
2. **Missing clinical CSV**: Pipeline expects CSV with same base filename as VCF.
3. **Sample ID mismatch**: Sample column name in VCF must exactly match `sample_id` in clinical CSV.
4. **Multi-sample VCFs**: Pipeline only processes single-sample VCFs by design.
5. **VRS API failures**: Network issues or invalid HGVSg may result in NA values for `vrs_id`. Pipeline continues processing.

## Authors

- Gabriele Bucci (bucci.gabriele@hsr.it) - ORCID: 0000-0001-9838-7204
- Michela Riba - ORCID: 0000-0003-4857-2027
- Muhammad Tajwar (contributor)
- Guido Scicolone (contributor)

## License

GPL-3
