# Platform-Specific Support Documentation

## Overview

This pipeline supports VCF files from multiple sequencing platforms with automatic detection and parameter optimization.

## Supported Platforms

### 1. Ion Torrent (Thermo Fisher Scientific)

**Detection Criteria:**
- Header contains: `##source="tvc X.XX-XX - Torrent Variant Caller"`
- Instrument field mentions: "IonTorrent", "S5", "Proton", "GeneStudio"

**Characteristics:**
- **Panel**: Oncomine Focus, Oncomine Comprehensive, etc.
- **Technology**: Semiconductor sequencing
- **Typical Coverage**: 500-5000x
- **Noise Profile**: Higher indel error rate, multi-allelic sites common
- **Specific Fields**: FAO (Filtered Allele Observation), FDP, FRO

**Optimized Parameters:**
```r
min_af = 0.01     # 1% - accounts for higher background noise
min_dp = 100      # Minimum depth
min_fao = 10      # Ion Torrent specific: minimum supporting reads
```

**Known Issues:**
- Multi-allelic variants with AF=0% for most alleles (artifacts)
- Solution: Automatically selects allele with highest FAO

### 2. Illumina

**Detection Criteria:**
- Header contains variant callers: `strelka`, `GATK`, `MuTect2`, `VarDict`, `VarScan`, etc.
- Platform field mentions: "Illumina", "NextSeq", "NovaSeq", "MiSeq"

**Characteristics:**
- **Panel**: Myriapod, TSO500, etc.
- **Technology**: Sequencing by synthesis
- **Typical Coverage**: 100-1000x
- **Noise Profile**: Lower error rate, cleaner multi-allelic calls
- **Specific Fields**: AD (Allelic Depth), DP

**Optimized Parameters:**
```r
min_af = 0.05     # 5% - stricter due to lower noise
min_dp = 50       # Lower depth acceptable
min_fao = NA      # Not used (FAO is Ion Torrent specific)
use_ad = TRUE     # Use AD field instead
```

**Planned Support** (not yet tested):
- Myriapod panel on Illumina
- TSO500 panel
- Custom panels

### 3. Future Platforms

**Under Consideration:**
- Oxford Nanopore (long-read)
- PacBio (long-read)
- BGI/MGI platforms

## Automatic Detection

The pipeline automatically detects the sequencing platform by analyzing VCF headers:

```r
# Automatic detection (default)
koios_vrs_pipeline(
    vcf_path = "sample.vcf",
    output_base_name = "sample_output",
    auto_detect_platform = TRUE  # default
)
```

### Detection Process

1. **Read VCF header** (first 200 lines)
2. **Check source field** for variant caller signatures
3. **Check platform/instrument fields**
4. **Identify panel** from header annotations
5. **Set optimal parameters** based on platform

### Detection Output

```
=== SEQUENCING PLATFORM DETECTION ===
Platform: IonTorrent
Variant Caller: TVC
Version: 5.12-27
Panel: Oncomine Focus

Recommended Parameters:
  min_af: 1.00%
  min_dp: 100
  min_fao: 10
```

## Manual Parameter Override

You can override auto-detected parameters:

```r
koios_vrs_pipeline(
    vcf_path = "sample.vcf",
    output_base_name = "sample_output",
    af_threshold = 0.02,  # Override to 2%
    min_dp = 200,         # Override to 200x
    min_fao = 20,         # Override to 20 reads
    auto_detect_platform = TRUE
)
```

## Disabling Auto-Detection

For manual control:

```r
koios_vrs_pipeline(
    vcf_path = "sample.vcf",
    output_base_name = "sample_output",
    af_threshold = 0.01,
    min_dp = 100,
    min_fao = 10,
    auto_detect_platform = FALSE
)
```

## Adding New Platform Support

To add support for a new platform:

### 1. Update `platform_detection.R`

Add detection logic:

```r
# --- YOUR PLATFORM DETECTION ---
your_platform_match <- grep("source.*your_caller", header_lines, value = TRUE)
if (length(your_platform_match) > 0) {
    platform_info$platform <- "YourPlatform"
    platform_info$caller <- "YourCaller"

    # Platform-specific parameters
    platform_info$params <- list(
        min_af = 0.XX,
        min_dp = XXX,
        min_fao = XX,
        use_fao = TRUE/FALSE,
        # Add custom parameters
    )

    return(platform_info)
}
```

### 2. Test with Real Data

```r
# Test detection
source("R/platform_detection.R")
platform_info <- detect_sequencing_platform("your_sample.vcf")
print_platform_summary(platform_info)

# Test full pipeline
result <- koios_vrs_pipeline(
    vcf_path = "your_sample.vcf",
    output_base_name = "test",
    auto_detect_platform = TRUE
)
```

### 3. Document Results

Update this file with:
- Detection criteria
- Platform characteristics
- Optimal parameters
- Known issues

## Real World Data (RWD) Examples

### Current RWD: Ion Torrent Oncomine Focus

**Location**: `RWD/Paz1_v1_Non-Filtered_2025-10-30_04_31_11.vcf`

**Platform Details:**
- Instrument: Ion S5
- Panel: Oncomine Focus
- Version: TVC 5.12-27
- Coverage: ~2000x median

**Processing Results:**
- Input: 543 variants
- After filtering: 24 high-quality SNP/INDEL variants
- Notable: EGFR L858R (p.Leu858Arg) correctly identified

### Future RWD: Illumina Myriapod

**Status**: Awaiting data

**Expected Characteristics:**
- Lower background noise
- Different INFO fields (AD instead of FAO)
- Cleaner multi-allelic representation

## Troubleshooting

### Platform Not Detected

**Symptom**: "Could not detect sequencing platform" warning

**Solutions:**
1. Check VCF header has `##source=` line
2. Manually specify parameters with `auto_detect_platform = FALSE`
3. Contact maintainer to add support for your platform

### Wrong Parameters Applied

**Symptom**: Too many/few variants pass filters

**Solutions:**
1. Verify detected platform with `detect_sequencing_platform()`
2. Override parameters manually
3. Check VCF header format

### Multi-allelic Issues

**Ion Torrent**: Multi-allelic sites with AF=0% are automatically filtered
**Illumina**: Multi-allelic sites usually represent real variants

If issues persist, adjust `min_fao` or `min_af` thresholds.

## Contact & Contributions

For questions or to add support for new platforms:
1. Open an issue with sample VCF header
2. Provide platform specifications
3. Share expected parameter ranges

## References

- Ion Torrent TVC: https://github.com/iontorrent/TS
- Illumina DRAGEN: https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html
- GATK Best Practices: https://gatk.broadinstitute.org/
