# KOIOS-VRS Pipeline with Platform Detection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17476991.svg)](https://doi.org/10.5281/zenodo.17476991)


Multi-platform genomic variant annotation pipeline integrating KOIOS, VRS, and Phenopackets v2.

## Features

- ✅ **Automatic Platform Detection**: Ion Torrent, Illumina, and more
- ✅ **Protein HGVS Annotations**: Both 3-letter and 1-letter formats
- ✅ **Multi-allelic Variant Filtering**: Automatic selection of highest-quality alleles
- ✅ **VRS Integration**: Variant Registration Service allele IDs
- ✅ **Phenopackets v2**: Standard genomic data exchange format
- ✅ **Interactive Clinical Data**: Prompts for missing metadata
- ✅ **hg19/hg38 Support**: Automatic liftover if needed

## Supported Platforms

| Platform | Status | Panel Examples | Auto-Detect |
|----------|--------|----------------|-------------|
| **Ion Torrent** | ✅ Tested | Oncomine Focus | ✅ Yes |
| **Illumina** | ⚠️ Ready | Myriapod, TSO500 | ✅ Yes |

See [PLATFORM_SUPPORT.md](PLATFORM_SUPPORT.md) for details.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Data Preparation](#input-data-preparation)
- [Usage](#usage)
- [Output Files](#output-files)
- [Advanced Options](#advanced-options)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

---

## Installation

### Requirements

- **R version**: ≥ 4.0.0
- **Operating System**: macOS, Linux, Windows

### Required R Packages

```r
# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "vcfR",
    "rtracklayer",
    "AnnotationHub",
    "GenomicRanges"
))

# Install from CRAN
install.packages(c(
    "httr",
    "jsonlite",
    "uuid"
))

# Install KOIOS from GitHub
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("CBIIT/KOIOS")
```

### Clone Repository

```bash
git clone https://github.com/YOUR_USERNAME/koios_vrs.git
cd koios_vrs
```

---

## Quick Start

### Basic Usage (Interactive R Session)

```r
source("R/koios_pipeline.R")

# Run pipeline with automatic platform detection
result <- koios_vrs_pipeline(
    vcf_path = "input/sample.vcf",
    output_base_name = "patient_001",
    ref_genome = "hg19"
)
```

### Command-Line Usage (Rscript)

```bash
Rscript -e "source('R/koios_pipeline.R'); \
  koios_vrs_pipeline(vcf_path='input/sample.vcf', \
                     output_base_name='patient_001', \
                     ref_genome='hg19')"
```

---

## Input Data Preparation

### 1. VCF File

Your VCF file should contain:
- Standard VCF format (v4.1 or later)
- Single sample (multi-sample VCFs not currently supported)
- Sequencing platform metadata in header (for auto-detection)

**Supported reference genomes**: `hg19` or `hg38` (automatic liftover if needed)

### 2. Clinical Metadata CSV

Create a CSV file with the same base name as your VCF:

**Location**: Same directory as VCF
**Naming**: `{vcf_basename}.csv`

**Example**: For `RWD/patient_001.vcf`, create `RWD/patient_001.csv`

**Required columns**:

```csv
sample_id,diagnosis_label,diagnosis_ncit_id,sex_label,sex_pato_id,age_of_onset_iso
patient_001,Lung Cancer,NCIT:C7377,MALE,PATO:0000384,P55Y
```

**Column descriptions**:
- `sample_id`: Patient/sample identifier (string)
- `diagnosis_label`: Human-readable diagnosis (e.g., "Lung Cancer")
- `diagnosis_ncit_id`: NCI Thesaurus ID (e.g., "NCIT:C7377")
- `sex_label`: "MALE" or "FEMALE"
- `sex_pato_id`: PATO ontology ID
  - MALE: `PATO:0000384`
  - FEMALE: `PATO:0000383`
- `age_of_onset_iso`: ISO 8601 duration format (e.g., "P55Y" = 55 years)

### Helper Script for CSV Creation

```r
# Create clinical CSV from command line
Rscript create_clinical_csv.R \
  input/patient_001.vcf \
  "Lung Cancer" \
  NCIT:C7377 \
  MALE \
  PATO:0000384 \
  55
```

Or run interactively:
```r
source("create_clinical_csv.R")
# Follow the prompts
```

---

## Usage

### Function Signature

```r
koios_vrs_pipeline(
    vcf_path,
    output_base_name,
    ref_genome = "hg38",
    af_threshold = NULL,
    min_dp = NULL,
    min_fao = NULL,
    auto_detect_platform = TRUE,
    input_dir = NULL,
    output_dir = "output"
)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vcf_path` | string | **required** | Path to input VCF file |
| `output_base_name` | string | **required** | Base name for output files |
| `ref_genome` | string | `"hg38"` | Reference genome: `"hg19"` or `"hg38"` |
| `af_threshold` | numeric | `NULL` | Minimum allele frequency (0-1). Auto-detected if NULL |
| `min_dp` | integer | `NULL` | Minimum depth of coverage. Auto-detected if NULL |
| `min_fao` | integer | `NULL` | Minimum filtered allele observations (Ion Torrent). Auto-detected if NULL |
| `auto_detect_platform` | logical | `TRUE` | Enable automatic platform detection |
| `input_dir` | string | `NULL` | Directory to search for clinical CSV |
| `output_dir` | string | `"output"` | Directory for output files |

### Auto-Detection Behavior

When `auto_detect_platform = TRUE` (default):

1. **Reads VCF header** to identify sequencing platform
2. **Detects variant caller** (TVC, GATK, Strelka, etc.)
3. **Sets optimal parameters** based on platform characteristics

**Ion Torrent defaults**:
- `af_threshold = 0.01` (1%)
- `min_dp = 100`
- `min_fao = 10`

**Illumina defaults**:
- `af_threshold = 0.05` (5%)
- `min_dp = 50`
- `min_fao = NA` (not used)

---

## Output Files

The pipeline generates 5-6 output files per run:

### 1. Phenopackets JSON (`*_Phenopackets.json`)

GA4GH Phenopackets v2 format with comprehensive variant annotations.

**Contents**:
- Patient metadata (diagnosis, sex, age)
- Variant interpretations with VRS IDs
- Genomic HGVS (e.g., `NC_000007.14:g.55191821_55191822delinsCG`)
- Protein HGVS 3-letter (e.g., `p.Leu858Arg`)
- Protein HGVS 1-letter (e.g., `p.L858R`)
- Allele frequency (e.g., `0.3302`)
- OMOP concept mappings

**Example structure**:
```json
{
  "id": "uuid",
  "subject": {
    "id": "patient_001",
    "sex": {"label": "MALE", "id": "PATO:0000384"}
  },
  "interpretations": [{
    "variantInterpretation": {
      "vrsAllele": {"id": "ga4gh:VA.xxx"},
      "extensions": [
        {"id": "hgvs_g", "value": "NC_000007.14:g.55191821_55191822delinsCG"},
        {"id": "allele_frequency", "value": 0.3302, "type": "number"},
        {"id": "hgvs_p_3letter", "value": "p.Leu858Arg"},
        {"id": "hgvs_p_1letter", "value": "p.L858R"}
      ]
    }
  }]
}
```

### 2. Note Variants CSV (`*_Note.csv`)

Clinically significant variants with OMOP concepts.

**Columns**:
- `CHROM`, `POS`, `REF`, `ALT`: Variant coordinates
- `hgvsg`: Genomic HGVS notation
- `hgvsp_3letter`: Protein HGVS 3-letter (e.g., p.Leu858Arg)
- `hgvsp_1letter`: Protein HGVS 1-letter (e.g., p.L858R)
- `allele_frequency`: Variant allele frequency (0-1)
- `concept_id`: OMOP concept ID
- `concept_name`: Human-readable concept name
- `vrs_id`: VRS allele identifier
- `URL`: VRS registry URL

### 3. VUS Variants CSV (`*_VUS.csv`)

Variants of Uncertain Significance (same format as Note.csv)

### 4. Gene Summary CSV (`*_GeneSummary.csv`)

Aggregated variant counts per gene.

**Columns**:
- `gene`: Gene symbol
- `concept_id`: OMOP concept ID
- `n_mutations`: Number of mutations in gene

### 5. Annotated VCF (`*_annotated.vcf`)

Original VCF with added annotations in INFO field:
- `VRS_ID`: VRS allele identifier
- `OMOP_CONCEPT`: OMOP concept ID
- `HGVSP_3L`: Protein HGVS 3-letter
- `HGVSP_1L`: Protein HGVS 1-letter

### 6. Intermediate Files (optional outputs)

- `*_hg38.vcf`: Lifted over VCF (if input was hg19)
- `*_hg38_preprocessed.vcf`: Filtered and preprocessed VCF
- `*_hg38_CNV.vcf`: Copy number variants (separated)

---

## Advanced Options

### Manual Parameter Override

Disable auto-detection and set custom filters:

```r
result <- koios_vrs_pipeline(
    vcf_path = "input/sample.vcf",
    output_base_name = "patient_strict",
    ref_genome = "hg19",
    auto_detect_platform = FALSE,
    af_threshold = 0.05,  # 5% minimum AF
    min_dp = 200,         # 200x minimum coverage
    min_fao = 20          # 20 filtered reads minimum
)
```

### Working with hg19 Data

The pipeline automatically lifts hg19 to hg38:

```r
result <- koios_vrs_pipeline(
    vcf_path = "input/legacy_sample.vcf",
    output_base_name = "legacy_patient",
    ref_genome = "hg19"  # Will be lifted to hg38 internally
)
```

### Specify Custom Directories

```r
result <- koios_vrs_pipeline(
    vcf_path = "data/vcf/sample.vcf",
    output_base_name = "analysis_001",
    input_dir = "data/clinical",  # Look for CSV here
    output_dir = "results/batch_01"
)
```

---

## Examples

### Example 1: Ion Torrent Oncomine Focus Panel

```r
source("R/koios_pipeline.R")

# VCF from Ion Torrent S5, Oncomine Focus panel
result <- koios_vrs_pipeline(
    vcf_path = "RWD/patient_lung_cancer.vcf",
    output_base_name = "LC_patient_001",
    ref_genome = "hg19"
)

# Auto-detects:
# - Platform: IonTorrent
# - Caller: TVC
# - Parameters: AF≥1%, DP≥100, FAO≥10
```

**Expected output**:
```
=== SEQUENCING PLATFORM DETECTION ===
Platform: IonTorrent
Variant Caller: TVC
Version: 5.12-27

Recommended Parameters:
  min_af: 1.00%
  min_dp: 100
  min_fao: 10
```

### Example 2: Illumina Panel (Future)

```r
# VCF from Illumina sequencer
result <- koios_vrs_pipeline(
    vcf_path = "illumina/sample.vcf",
    output_base_name = "ILL_patient_001",
    ref_genome = "hg38"
)

# Auto-detects:
# - Platform: Illumina
# - Caller: GATK/Strelka/etc
# - Parameters: AF≥5%, DP≥50
```

### Example 3: Batch Processing

```r
# Process multiple VCFs
vcf_files <- list.files("input", pattern = "\\.vcf$", full.names = TRUE)

for (vcf_file in vcf_files) {
    base_name <- tools::file_path_sans_ext(basename(vcf_file))

    result <- koios_vrs_pipeline(
        vcf_path = vcf_file,
        output_base_name = base_name,
        ref_genome = "hg19",
        output_dir = "batch_results"
    )

    cat(sprintf("Processed: %s\n", base_name))
}
```

---

## Troubleshooting

### Issue 1: Clinical CSV Not Found

**Error**:
```
ERROR: Clinical metadata CSV file not found.
```

**Solution**:
1. Create CSV with same base name as VCF
2. Place in same directory as VCF
3. Or use `create_clinical_csv.R` helper script

```bash
Rscript create_clinical_csv.R input/sample.vcf "Lung Cancer" NCIT:C7377 MALE PATO:0000384 55
```

### Issue 2: Platform Not Detected

**Warning**:
```
Could not detect sequencing platform. Using conservative default parameters.
```

**Solution**:
- Check VCF header contains platform information
- Manually specify parameters with `auto_detect_platform = FALSE`
- See [PLATFORM_SUPPORT.md](PLATFORM_SUPPORT.md) for detection criteria

### Issue 3: No Variants Pass Filters

**Warning**:
```
No variants remained after filtering.
```

**Possible causes**:
- AF threshold too high
- DP threshold too high
- FAO threshold too high (Ion Torrent)

**Solution**: Lower thresholds manually:
```r
result <- koios_vrs_pipeline(
    vcf_path = "input/sample.vcf",
    output_base_name = "lenient",
    af_threshold = 0.005,  # 0.5%
    min_dp = 50,
    min_fao = 5
)
```

### Issue 4: VRS API Timeout

**Error**:
```
Failed to fetch VRS ID for variant
```

**Cause**: VRS registry temporarily unavailable or slow

**Solution**:
- Check internet connection
- Retry pipeline (results are cached)
- VRS IDs will be `NA` for failed queries

### Issue 5: Multi-allelic Variants

**Expected behavior**: Pipeline automatically selects highest-quality allele (highest FAO/AD)

If you see unexpected duplicates:
- Check `*_preprocessed.vcf` for multi-allelic sites
- Verify FAO/AD values in original VCF
- Report issue with example VCF

---

## Citation

If you use this pipeline in your research, please cite:

```
@software{koios_vrs_pipeline,
  author = {Your Name},
  title = {KOIOS-VRS Pipeline: Multi-platform Genomic Variant Annotation},
  year = {2024},
  doi = {10.5281/zenodo.17476991},
  url = {https://github.com/YOUR_USERNAME/koios_vrs}
}
```

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

## License

This project is licensed under the MIT License.

---

## Support

For issues, questions, or feature requests:
- Open an issue on GitHub
- See [PLATFORM_SUPPORT.md](PLATFORM_SUPPORT.md) for platform-specific guidance

---

## Acknowledgments

- **KOIOS**: OMOP-based genomic variant annotation
- **VRS**: GA4GH Variant Representation Specification
- **Phenopackets**: GA4GH standard for phenotypic data
