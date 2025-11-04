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

## Quick Start

```r
source("R/koios_pipeline.R")

# Auto-detection (recommended)
result <- koios_vrs_pipeline(
    vcf_path = "input/sample.vcf",
    output_base_name = "patient_A",
    ref_genome = "hg19"
)
```

See full documentation in this file.
