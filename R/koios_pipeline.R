# koios_pipeline.R

#' Preprocess VCF: Normalize, Filter, and Separate Variants
#'
#' Preprocesses VCF to optimize for KOIOS processing by:
#' - Filtering low quality variants (AF < 1%, DP < 100)
#' - Removing undefined genotypes (./.)
#' - Selecting most representative allele for multi-allelic variants
#' - Separating SNP/INDEL from CNV variants
#'
#' @param vcf_obj vcfR object to preprocess
#' @param min_af Minimum allele frequency threshold (default: 0.01)
#' @param min_dp Minimum read depth threshold (default: 100)
#' @param min_fao Minimum filtered allele observation (default: 10)
#' @return List with 'snp_indel' vcfR object and 'cnv' vcfR object
preprocess_vcf <- function(vcf_obj, min_af = 0.01, min_dp = 100, min_fao = 10) {

    cat("\n=== VCF PREPROCESSING ===\n")

    # Get data
    fix <- vcf_obj@fix
    gt <- vcf_obj@gt
    n_initial <- nrow(fix)
    cat(sprintf("Initial variants: %d\n", n_initial))

    # 1. Filter undefined genotypes (./.)
    cat("Step 1: Filtering undefined genotypes (./.)\n")
    sample_gt <- gt[, 2]  # Assuming single sample VCF
    has_valid_gt <- !grepl("^\\./\\.$", sample_gt)
    keep_idx <- has_valid_gt
    cat(sprintf("  Removed %d variants with undefined genotypes\n", sum(!has_valid_gt)))

    # 2. Extract and filter by AF and DP
    cat("Step 2: Filtering by AF and DP\n")
    info <- fix[keep_idx, "INFO"]

    # Extract AF values
    af_vals <- sapply(info, function(x) {
        af_match <- regmatches(x, regexec("AF=([0-9.,]+)", x))
        if (length(af_match[[1]]) > 1) {
            # Get all AF values (comma-separated for multi-allelic)
            afs <- as.numeric(strsplit(af_match[[1]][2], ",")[[1]])
            return(max(afs, na.rm = TRUE))  # Take maximum AF
        }
        return(NA_real_)
    })

    # Extract DP values
    dp_vals <- sapply(info, function(x) {
        dp_match <- regmatches(x, regexec("DP=([0-9]+)", x))
        if (length(dp_match[[1]]) > 1) {
            return(as.numeric(dp_match[[1]][2]))
        }
        return(NA_real_)
    })

    # Extract FAO values (maximum across alleles for multi-allelic)
    fao_vals <- sapply(info, function(x) {
        fao_match <- regmatches(x, regexec("FAO=([0-9,]+)", x))
        if (length(fao_match[[1]]) > 1) {
            faos <- as.numeric(strsplit(fao_match[[1]][2], ",")[[1]])
            return(max(faos, na.rm = TRUE))
        }
        return(NA_real_)
    })

    # Apply filters
    pass_af <- is.na(af_vals) | af_vals >= min_af
    pass_dp <- is.na(dp_vals) | dp_vals >= min_dp
    pass_fao <- is.na(fao_vals) | fao_vals >= min_fao
    filter_pass <- pass_af & pass_dp & pass_fao

    cat(sprintf("  Removed %d variants with AF < %.1f%%\n", sum(!pass_af), min_af * 100))
    cat(sprintf("  Removed %d variants with DP < %d\n", sum(!pass_dp), min_dp))
    cat(sprintf("  Removed %d variants with FAO < %d\n", sum(!pass_fao), min_fao))

    keep_idx[keep_idx] <- filter_pass

    # 3. Separate CNV from SNP/INDEL
    cat("Step 3: Separating CNV from SNP/INDEL\n")
    alt_vals <- fix[keep_idx, "ALT"]
    is_cnv <- grepl("<CNV>", alt_vals, fixed = TRUE)

    cat(sprintf("  Found %d CNV variants\n", sum(is_cnv)))
    cat(sprintf("  Found %d SNP/INDEL variants\n", sum(!is_cnv)))

    # 4. Handle multi-allelic variants (for SNP/INDEL only)
    cat("Step 4: Selecting most representative allele for multi-allelic variants\n")

    snp_indel_idx <- which(keep_idx)[!is_cnv]
    cnv_idx <- which(keep_idx)[is_cnv]

    if (length(snp_indel_idx) > 0) {
        # For SNP/INDEL: select allele with highest FAO or AF and remove other alleles
        snp_fix <- fix[snp_indel_idx, , drop = FALSE]
        snp_gt <- gt[snp_indel_idx, , drop = FALSE]

        # Process each variant to keep only the best allele
        n_multi_allelic <- 0
        for (i in seq_len(nrow(snp_fix))) {
            alt_str <- snp_fix[i, "ALT"]
            alts <- strsplit(alt_str, ",")[[1]]

            # If multi-allelic, select best allele
            if (length(alts) > 1) {
                n_multi_allelic <- n_multi_allelic + 1
                info_str <- snp_fix[i, "INFO"]

                # Extract FAO values for each allele
                fao_match <- regmatches(info_str, regexec("FAO=([0-9,]+)", info_str))
                if (length(fao_match[[1]]) > 1) {
                    faos <- as.numeric(strsplit(fao_match[[1]][2], ",")[[1]])

                    # Find allele with maximum FAO
                    if (length(faos) == length(alts)) {
                        best_idx <- which.max(faos)
                        best_alt <- alts[best_idx]

                        # Update ALT field to contain only the best allele
                        snp_fix[i, "ALT"] <- best_alt

                        # Note: We keep the INFO field as-is because it's complex to modify
                        # KOIOS will handle the multi-allelic INFO fields
                    }
                }
            }
        }

        # Now remove duplicate positions (keeping first which has highest FAO due to sorting)
        # Extract FAO for sorting
        fao_vals <- sapply(seq_len(nrow(snp_fix)), function(i) {
            info_str <- snp_fix[i, "INFO"]
            fao_match <- regmatches(info_str, regexec("FAO=([0-9,]+)", info_str))
            if (length(fao_match[[1]]) > 1) {
                faos <- as.numeric(strsplit(fao_match[[1]][2], ",")[[1]])
                return(max(faos, na.rm = TRUE))
            }
            return(0)
        })

        # Group by position and select best allele
        pos_key <- paste(snp_fix[, "CHROM"], snp_fix[, "POS"])
        n_before_dedup <- length(snp_indel_idx)

        # Order by FAO descending
        order_idx <- order(pos_key, -fao_vals)
        snp_fix <- snp_fix[order_idx, , drop = FALSE]
        snp_gt <- snp_gt[order_idx, , drop = FALSE]
        pos_key <- pos_key[order_idx]

        # Keep first (best) of each position
        keep_dedup <- !duplicated(pos_key)
        snp_fix <- snp_fix[keep_dedup, , drop = FALSE]
        snp_gt <- snp_gt[keep_dedup, , drop = FALSE]

        n_after_dedup <- nrow(snp_fix)
        cat(sprintf("  Processed %d multi-allelic variants\n", n_multi_allelic))
        cat(sprintf("  Removed %d duplicate positions (kept allele with highest FAO)\n",
                    n_before_dedup - n_after_dedup))

        # Create SNP/INDEL vcfR object
        snp_vcf <- vcf_obj
        snp_vcf@fix <- snp_fix
        snp_vcf@gt <- snp_gt
    } else {
        snp_vcf <- NULL
    }

    # Create CNV vcfR object
    if (length(cnv_idx) > 0) {
        cnv_vcf <- vcf_obj
        cnv_vcf@fix <- fix[cnv_idx, , drop = FALSE]
        cnv_vcf@gt <- gt[cnv_idx, , drop = FALSE]
    } else {
        cnv_vcf <- NULL
    }

    cat(sprintf("\nPreprocessing complete:\n"))
    cat(sprintf("  Initial: %d variants\n", n_initial))
    cat(sprintf("  Final SNP/INDEL: %d variants\n", if(!is.null(snp_vcf)) nrow(snp_vcf@fix) else 0))
    cat(sprintf("  Final CNV: %d variants\n", if(!is.null(cnv_vcf)) nrow(cnv_vcf@fix) else 0))
    cat(sprintf("  Total filtered: %d variants\n\n",
                n_initial - ifelse(!is.null(snp_vcf), nrow(snp_vcf@fix), 0) -
                           ifelse(!is.null(cnv_vcf), nrow(cnv_vcf@fix), 0)))

    return(list(snp_indel = snp_vcf, cnv = cnv_vcf))
}

#' Convert Three-Letter to Single-Letter Amino Acid Code
#'
#' Converts three-letter amino acid codes to single-letter codes
#'
#' @param three_letter Three-letter amino acid code
#' @return Single-letter amino acid code
three_to_one_aa <- function(three_letter) {
    aa_map <- c(
        "Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D",
        "Cys" = "C", "Gln" = "Q", "Glu" = "E", "Gly" = "G",
        "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K",
        "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S",
        "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V",
        "Ter" = "*", "Stop" = "*"
    )

    if (three_letter %in% names(aa_map)) {
        return(aa_map[three_letter])
    }
    return(NA_character_)
}

#' Convert Single-Letter to Three-Letter Amino Acid Code
#'
#' Converts single-letter amino acid codes to three-letter codes
#'
#' @param one_letter Single-letter amino acid code
#' @return Three-letter amino acid code
one_to_three_aa <- function(one_letter) {
    aa_map <- c(
        "A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp",
        "C" = "Cys", "Q" = "Gln", "E" = "Glu", "G" = "Gly",
        "H" = "His", "I" = "Ile", "L" = "Leu", "K" = "Lys",
        "M" = "Met", "F" = "Phe", "P" = "Pro", "S" = "Ser",
        "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val",
        "*" = "Ter"
    )

    if (one_letter %in% names(aa_map)) {
        return(aa_map[one_letter])
    }
    return(NA_character_)
}

#' Extract Protein HGVS Annotations from VCF FUNC Field
#'
#' Extracts protein annotations from the FUNC field in VCF INFO column.
#' Generates both three-letter and single-letter HGVS protein notations.
#'
#' @param vcf_df Data frame with VCF data including INFO column
#' @return Data frame with added hgvsp_3letter and hgvsp_1letter columns
extract_protein_hgvs <- function(vcf_df) {
    if (!"INFO" %in% names(vcf_df)) {
        vcf_df$hgvsp_3letter <- NA_character_
        vcf_df$hgvsp_1letter <- NA_character_
        vcf_df$allele_frequency <- NA_real_
        return(vcf_df)
    }

    hgvsp_3letter <- character(nrow(vcf_df))
    hgvsp_1letter <- character(nrow(vcf_df))
    allele_frequency <- numeric(nrow(vcf_df))

    for (i in seq_len(nrow(vcf_df))) {
        info_field <- vcf_df$INFO[i]

        if (is.na(info_field) || !nzchar(info_field)) {
            hgvsp_3letter[i] <- NA_character_
            hgvsp_1letter[i] <- NA_character_
            allele_frequency[i] <- NA_real_
            next
        }

        # Extract AF (Allele Frequency)
        af_match <- regmatches(info_field, regexec("AF=([0-9.,]+)", info_field))
        if (length(af_match[[1]]) > 1) {
            # Get all AF values (comma-separated for multi-allelic)
            afs <- as.numeric(strsplit(af_match[[1]][2], ",")[[1]])
            allele_frequency[i] <- max(afs, na.rm = TRUE)  # Take maximum AF
        } else {
            allele_frequency[i] <- NA_real_
        }

        # Extract FUNC field which contains protein annotation
        # Use a greedy pattern to capture everything between FUNC=[ and the closing ]
        func_match <- regmatches(info_field, regexpr("FUNC=\\[.*\\]", info_field))

        if (length(func_match) > 0) {
            func_str <- func_match[1]

            # Extract protein field (e.g., 'protein':'p.Leu858Arg')
            protein_match <- regmatches(func_str, regexpr("'protein':'(p\\.[^']+)'", func_str, perl = TRUE))

            if (length(protein_match) > 0) {
                # Extract p.Leu858Arg format (remove 'protein':' prefix and trailing ')
                protein_str <- sub("'protein':'", "", protein_match[1])
                protein_str <- sub("'$", "", protein_str)  # Remove trailing quote if present

                # This is already in three-letter format
                hgvsp_3letter[i] <- protein_str

                # Convert to single-letter format
                # Parse the protein change (e.g., p.Leu858Arg)
                if (grepl("^p\\.", protein_str)) {
                    # Extract components using regex
                    # Pattern: p.Xxx###Yyy or p.Xxx###fs, p.Xxx###del, etc.

                    # Handle substitutions (e.g., p.Leu858Arg)
                    if (grepl("^p\\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}$", protein_str)) {
                        ref_aa <- sub("^p\\.([A-Z][a-z]{2}).*", "\\1", protein_str)
                        pos <- sub("^p\\.[A-Z][a-z]{2}([0-9]+).*", "\\1", protein_str)
                        alt_aa <- sub("^p\\.[A-Z][a-z]{2}[0-9]+(.*)", "\\1", protein_str)

                        ref_aa_1 <- three_to_one_aa(ref_aa)
                        alt_aa_1 <- three_to_one_aa(alt_aa)

                        if (!is.na(ref_aa_1) && !is.na(alt_aa_1)) {
                            hgvsp_1letter[i] <- sprintf("p.%s%s%s", ref_aa_1, pos, alt_aa_1)
                        } else {
                            hgvsp_1letter[i] <- NA_character_
                        }
                    } else {
                        # Keep as is for complex variants (del, ins, dup, fs, etc.)
                        hgvsp_1letter[i] <- protein_str
                    }
                } else {
                    hgvsp_1letter[i] <- NA_character_
                }
            } else {
                hgvsp_3letter[i] <- NA_character_
                hgvsp_1letter[i] <- NA_character_
            }
        } else {
            hgvsp_3letter[i] <- NA_character_
            hgvsp_1letter[i] <- NA_character_
        }
    }

    vcf_df$hgvsp_3letter <- hgvsp_3letter
    vcf_df$hgvsp_1letter <- hgvsp_1letter
    vcf_df$allele_frequency <- allele_frequency

    return(vcf_df)
}

#' Prompt User for Clinical Metadata Interactively
#'
#' Prompts the user to enter clinical metadata for a sample when CSV file is missing.
#' Works in both interactive R sessions and Rscript mode.
#'
#' @param sample_id The sample ID from the VCF file
#' @return A data frame with clinical metadata
prompt_clinical_data <- function(sample_id) {
    cat("\n=== CLINICAL METADATA REQUIRED ===\n")
    cat(sprintf("Sample ID: %s\n\n", sample_id))

    # Check if running interactively
    is_interactive <- interactive()

    if (!is_interactive) {
        cat("ERROR: Clinical metadata CSV file not found.\n")
        cat(sprintf("Please create a CSV file with the following format:\n"))
        cat(sprintf("Location: Same directory as VCF, named '{vcf_basename}.csv'\n\n"))
        cat("Required columns:\n")
        cat("  sample_id,diagnosis_label,diagnosis_ncit_id,sex_label,sex_pato_id,age_of_onset_iso\n\n")
        cat("Example:\n")
        cat(sprintf("  %s,Lung Cancer,NCIT:C7377,MALE,PATO:0000384,P45Y\n\n", sample_id))
        cat("Common sex values:\n")
        cat("  MALE: PATO:0000384\n")
        cat("  FEMALE: PATO:0000383\n\n")
        cat("Age format: P{years}Y (e.g., P45Y for 45 years old)\n\n")
        cat("Once the CSV file is created, run the pipeline again.\n")
        stop("Clinical metadata CSV file required but not found. Cannot prompt in non-interactive mode.")
    }

    # Interactive mode - use readline
    cat("CSV file not found. Please enter clinical data interactively:\n\n")

    # Diagnosis
    diagnosis_label <- readline(prompt = "Enter diagnosis label (e.g., 'Lung Cancer'): ")
    if (!nzchar(diagnosis_label)) {
        stop("Diagnosis label is required")
    }

    diagnosis_ncit_id <- readline(prompt = "Enter diagnosis NCIT ID (e.g., 'NCIT:C7377'): ")
    if (!nzchar(diagnosis_ncit_id)) {
        stop("Diagnosis NCIT ID is required")
    }

    # Sex
    cat("\nSex options:\n")
    cat("  1. MALE (PATO:0000384)\n")
    cat("  2. FEMALE (PATO:0000383)\n")
    sex_choice <- readline(prompt = "Select sex (1 or 2): ")

    if (sex_choice == "1") {
        sex_label <- "MALE"
        sex_pato_id <- "PATO:0000384"
    } else if (sex_choice == "2") {
        sex_label <- "FEMALE"
        sex_pato_id <- "PATO:0000383"
    } else {
        sex_label <- readline(prompt = "Enter sex label: ")
        sex_pato_id <- readline(prompt = "Enter sex PATO ID: ")
    }

    # Age of onset
    age_of_onset <- readline(prompt = "Enter age of onset in years (e.g., '45'): ")
    if (!nzchar(age_of_onset)) {
        stop("Age of onset is required")
    }
    age_of_onset_iso <- sprintf("P%sY", age_of_onset)

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

    cat("\n=== CLINICAL METADATA SUMMARY ===\n")
    print(clinical_data)
    cat("\n")

    return(clinical_data)
}

#' KOIOS and VRS Pipeline for Cancer Variants
#'
#' Processes a VCF file, enriches it with KOIOS (HGVSg, ClinGen, OMOP),
#' queries the VRS API for VRS IDs, filters variants based on frequency,
#' and generates Phenopackets v2 JSON output along with VUS/Note CSV files.
#'
#' @param vcf_path Path to the input VCF file (hg19 or hg38). Can be absolute or relative to input_dir.
#' @param output_base_name Base name for output files (e.g., 'patient_A' produces 'patient_A_VUS.csv', etc.).
#' @param ref_genome Reference genome of input VCF: "hg19" or "hg38" (default: "hg38").
#' @param af_threshold Allele Frequency threshold for filtering (default: 0.01 or 1%). Variants with population AF > threshold are excluded.
#' @param min_dp Minimum read depth for variant quality filtering (default: 100).
#' @param min_fao Minimum filtered allele observation for variant quality filtering (default: 10).
#' @param input_dir Input directory containing VCF and clinical CSV files (default: "input").
#' @param output_dir Output directory for all generated files (default: "output").
#' @param default_vus_id OMOP Concept ID for Variant of Unknown Significance (default: 1028197L).
#' @return A list containing the VUS and Note dataframes, invisibly. Writes output files to disk.
#' @export
koios_vrs_pipeline <- function(vcf_path, output_base_name, ref_genome = "hg38", af_threshold = NULL, min_dp = NULL, min_fao = NULL,
                                input_dir = "input", output_dir = "output", default_vus_id = 1028197L, auto_detect_platform = TRUE) {

    # --- 1. SETUP AND INPUT VALIDATION ---

    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat(sprintf("Created output directory: %s\n", output_dir))
    }

    # Validate reference genome parameter
    if (!ref_genome %in% c("hg19", "hg38")) {
        stop("ref_genome must be either 'hg19' or 'hg38'. Got: ", ref_genome)
    }

    # Resolve input file path
    if (!file.exists(vcf_path)) {
        vcf_path_alt <- file.path(input_dir, vcf_path)
        if (file.exists(vcf_path_alt)) {
            vcf_path <- vcf_path_alt
        } else {
            stop(sprintf("VCF file not found: %s (also checked %s)", vcf_path, vcf_path_alt))
        }
    }

    # --- 1.1. PLATFORM DETECTION ---
    if (auto_detect_platform) {
        source("R/platform_detection.R")
        platform_info <- detect_sequencing_platform(vcf_path)
        print_platform_summary(platform_info)

        # Use platform-specific parameters if not explicitly provided
        if (is.null(af_threshold)) {
            af_threshold <- platform_info$params$min_af
        }
        if (is.null(min_dp)) {
            min_dp <- platform_info$params$min_dp
        }
        if (is.null(min_fao)) {
            min_fao <- ifelse(is.na(platform_info$params$min_fao), 10, platform_info$params$min_fao)
        }
    } else {
        # Use defaults if auto-detection is disabled
        if (is.null(af_threshold)) af_threshold <- 0.01
        if (is.null(min_dp)) min_dp <- 100
        if (is.null(min_fao)) min_fao <- 10
    }

    # 1.2. Get Sample ID from VCF (must be done BEFORE clinical CSV check)
    vcf_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
    sample_names <- colnames(vcf_obj@gt)[-1]
    if (length(sample_names) != 1) {
        stop("VCF must contain exactly one sample column. Found: ", length(sample_names))
    }
    sample_id <- sample_names[1]

    # 1.3. Find and load clinical CSV (same base name as VCF)
    vcf_basename <- basename(vcf_path)
    csv_basename <- sub("\\.vcf$", ".csv", vcf_basename, ignore.case = TRUE)
    clinical_csv_path <- file.path(dirname(vcf_path), csv_basename)

    if (!file.exists(clinical_csv_path)) {
        cat(sprintf("Clinical metadata file not found: %s\n", clinical_csv_path))
        clinical_data <- prompt_clinical_data(sample_id)
        # Write to CSV for future use
        utils::write.csv(clinical_data, file = clinical_csv_path, row.names = FALSE, quote = TRUE)
        cat(sprintf("Clinical metadata saved to: %s\n", clinical_csv_path))
    } else {
        # Load Clinical Metadata
        clinical_data <- utils::read.csv(clinical_csv_path, stringsAsFactors = FALSE)
        clinical_data <- clinical_data[clinical_data$sample_id == sample_id, ]

        if (nrow(clinical_data) != 1) {
            stop("Clinical metadata not found or does not correspond to the VCF sample ID in the CSV.")
        }
    }

    # --- 1.4. LIFTOVER FROM HG19 TO HG38 IF NEEDED ---

    if (ref_genome == "hg19") {
        cat("Input VCF is hg19. Performing liftOver to hg38...\n")

        # Create temporary hg38 VCF file path in output directory
        vcf_hg38_filename <- sub("\\.vcf$", "_hg38.vcf", basename(vcf_path), ignore.case = TRUE)
        vcf_hg38_path <- file.path(output_dir, vcf_hg38_filename)

        # Load required libraries for liftover
        if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("Package 'rtracklayer' is required for hg19 to hg38 liftover. Install with: BiocManager::install('rtracklayer')")
        }
        if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
            stop("Package 'GenomicRanges' is required for hg19 to hg38 liftover. Install with: BiocManager::install('GenomicRanges')")
        }
        if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
            stop("Package 'AnnotationHub' is required for liftover chain file. Install with: BiocManager::install('AnnotationHub')")
        }

        # Get liftover chain file from AnnotationHub
        ah <- AnnotationHub::AnnotationHub()
        chain <- AnnotationHub::query(ah, c("hg19", "hg38", "chain"))
        chain_file <- chain[["AH14150"]]  # hg19ToHg38.over.chain.gz

        # Read VCF and convert to GRanges
        vcf_data <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
        vcf_fix <- vcf_data@fix  # Access fix slot directly for all columns including INFO
        vcf_gt <- vcf_data@gt

        # Create GRanges object for liftover
        # Ensure chromosome names have 'chr' prefix (add if missing)
        chrom_names <- vcf_fix[, "CHROM"]
        if (!grepl("^chr", chrom_names[1])) {
            chrom_names <- paste0("chr", chrom_names)
        }

        gr_hg19 <- GenomicRanges::GRanges(
            seqnames = chrom_names,
            ranges = IRanges::IRanges(
                start = as.numeric(vcf_fix[, "POS"]),
                end = as.numeric(vcf_fix[, "POS"])
            )
        )

        # Add metadata (store original index to reconstruct VCF later)
        GenomicRanges::mcols(gr_hg19) <- data.frame(
            original_index = seq_len(nrow(vcf_fix)),
            stringsAsFactors = FALSE
        )

        # Perform liftover
        gr_hg38_list <- rtracklayer::liftOver(gr_hg19, chain_file)

        # Extract successfully lifted variants (keep only one-to-one mappings)
        lifted_lengths <- lengths(gr_hg38_list)
        successfully_lifted <- which(lifted_lengths == 1)

        if (length(successfully_lifted) == 0) {
            stop("No variants could be successfully lifted from hg19 to hg38.")
        }

        cat(sprintf("Successfully lifted %d/%d variants from hg19 to hg38.\n",
                    length(successfully_lifted), length(gr_hg19)))

        # Extract lifted coordinates
        gr_hg38 <- unlist(gr_hg38_list[successfully_lifted])

        # Reconstruct VCF in hg38
        new_chrom <- as.character(GenomicRanges::seqnames(gr_hg38))
        new_chrom <- sub("^chr", "", new_chrom)  # Remove chr prefix for consistency
        new_pos <- GenomicRanges::start(gr_hg38)

        # Create new VCF fix data with updated coordinates
        lifted_indices <- GenomicRanges::mcols(gr_hg38)$original_index
        lifted_vcf_fix <- vcf_fix[lifted_indices, , drop = FALSE]
        lifted_vcf_fix[, "CHROM"] <- new_chrom
        lifted_vcf_fix[, "POS"] <- as.character(new_pos)

        # Create new VCF gt data (subset to lifted variants)
        lifted_vcf_gt <- vcf_gt[lifted_indices, , drop = FALSE]

        # Write new hg38 VCF
        vcf_hg38_obj <- vcf_data
        vcf_hg38_obj@fix <- lifted_vcf_fix
        vcf_hg38_obj@gt <- lifted_vcf_gt

        # Write uncompressed VCF
        temp_hg38_gz_path <- paste0(vcf_hg38_path, ".gz")
        vcfR::write.vcf(vcf_hg38_obj, file = temp_hg38_gz_path)
        system(sprintf("gunzip -f '%s'", temp_hg38_gz_path))
        cat(sprintf("Lifted VCF written to: %s\n", vcf_hg38_path))

        # Use the lifted VCF for the rest of the pipeline
        vcf_path <- vcf_hg38_path
        # Reload vcf_obj with lifted coordinates
        vcf_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
    }

    # --- 1.5. VCF PREPROCESSING ---

    # Preprocess VCF: filter quality, remove multi-allelic, separate CNV
    preprocessed <- preprocess_vcf(vcf_obj, min_af = af_threshold, min_dp = min_dp, min_fao = min_fao)

    if (is.null(preprocessed$snp_indel) || nrow(preprocessed$snp_indel@fix) == 0) {
        warning("No SNP/INDEL variants remained after preprocessing. Skipping KOIOS processing.")
        return(invisible(list(vus_df = data.frame(), note_df = data.frame())))
    }

    # Write preprocessed VCF for KOIOS in output directory (uncompressed)
    preprocessed_vcf_filename <- sub("\\.vcf$", "_preprocessed.vcf", basename(vcf_path), ignore.case = TRUE)
    preprocessed_vcf_path <- file.path(output_dir, preprocessed_vcf_filename)

    # Write uncompressed VCF by using a temp file and then gunzip
    temp_gz_path <- paste0(preprocessed_vcf_path, ".gz")
    vcfR::write.vcf(preprocessed$snp_indel, file = temp_gz_path)
    system(sprintf("gunzip -f '%s'", temp_gz_path))
    cat(sprintf("Preprocessed VCF written to: %s\n\n", preprocessed_vcf_path))

    # Optional: Write CNV to separate file for future processing (uncompressed)
    if (!is.null(preprocessed$cnv) && nrow(preprocessed$cnv@fix) > 0) {
        cnv_vcf_filename <- sub("\\.vcf$", "_CNV.vcf", basename(vcf_path), ignore.case = TRUE)
        cnv_vcf_path <- file.path(output_dir, cnv_vcf_filename)
        temp_cnv_gz_path <- paste0(cnv_vcf_path, ".gz")
        vcfR::write.vcf(preprocessed$cnv, file = temp_cnv_gz_path)
        system(sprintf("gunzip -f '%s'", temp_cnv_gz_path))
        cat(sprintf("CNV variants written to: %s\n\n", cnv_vcf_path))
    }

    # --- 2. EXTRACT PROTEIN HGVS ANNOTATIONS FROM ORIGINAL VCF ---

    # Extract protein annotations from the preprocessed VCF BEFORE KOIOS processing
    # (KOIOS processVCF does not preserve the INFO field)
    cat("\n=== EXTRACTING PROTEIN HGVS ANNOTATIONS ===\n")
    preprocessed_vcf_obj <- vcfR::read.vcfR(preprocessed_vcf_path, verbose = FALSE)

    # Create temporary data frame with VCF data for annotation extraction
    # Note: After preprocessing, ALT field should contain only one allele per variant
    # but we still check for multi-allelic variants just in case
    temp_vcf_df <- data.frame(
        CHROM = preprocessed_vcf_obj@fix[, "CHROM"],
        POS = preprocessed_vcf_obj@fix[, "POS"],
        REF = preprocessed_vcf_obj@fix[, "REF"],
        ALT = preprocessed_vcf_obj@fix[, "ALT"],
        INFO = preprocessed_vcf_obj@fix[, "INFO"],
        stringsAsFactors = FALSE
    )

    # Extract protein annotations
    temp_vcf_df <- extract_protein_hgvs(temp_vcf_df)
    n_with_protein <- sum(!is.na(temp_vcf_df$hgvsp_3letter))
    cat(sprintf("Found protein annotations for %d/%d variants (post-explosion)\n\n", n_with_protein, nrow(temp_vcf_df)))

    # --- 3. KOIOS PROCESSING ---

    cat("=== KOIOS PROCESSING ===\n")
    concepts <- KOIOS::loadConcepts()
    vcf <- KOIOS::loadVCF(userVCF = preprocessed_vcf_path)  # Use preprocessed VCF
    ref.df <- KOIOS::loadReference("hg38")

    vcf.df <- KOIOS::processVCF(vcf)
    vcf.df <- KOIOS::generateHGVSG(vcf = vcf.df, ref = ref.df)
    vcf.df <- KOIOS::processClinGen(vcf.df, ref = "hg38", progressBar = FALSE)
    vcf.df <- KOIOS::addConcepts(vcf.df, concepts, returnAll = FALSE)

    default_vus_name <- "Variant of unknown significance"
    vcf.df$concept_id[is.na(vcf.df$concept_id)] <- default_vus_id
    vcf.df$concept_name[is.na(vcf.df$concept_name) | !nzchar(vcf.df$concept_name)] <- default_vus_name

    # --- 3.1. MERGE PROTEIN ANNOTATIONS WITH KOIOS DATA ---
    cat("\n=== MERGING PROTEIN ANNOTATIONS ===\n")

    # Merge by chromosome, position, ref, and alt
    vcf.df <- merge(
        vcf.df,
        temp_vcf_df[, c("CHROM", "POS", "REF", "ALT", "hgvsp_3letter", "hgvsp_1letter", "allele_frequency")],
        by = c("CHROM", "POS", "REF", "ALT"),
        all.x = TRUE,
        sort = FALSE
    )

    n_merged <- sum(!is.na(vcf.df$hgvsp_3letter))
    cat(sprintf("Merged protein annotations for %d/%d variants\n\n", n_merged, nrow(vcf.df)))

    # --- 4. POST-KOIOS FILTERING ---
    # Note: Quality filtering, CNV removal, and multi-allelic resolution already done in preprocessing

    # Optional: Additional population frequency filtering based on gnomAD
    freq_col <- "koios_gnomAD_AF"
    if (freq_col %in% names(vcf.df)) {
        vcf.df[[freq_col]] <- suppressWarnings(as.numeric(vcf.df[[freq_col]]))
        n_before_gnomad <- nrow(vcf.df)
        exclude_freq <- !is.na(vcf.df[[freq_col]]) & vcf.df[[freq_col]] > af_threshold
        variants_for_vrs <- vcf.df[!exclude_freq, ]
        cat(sprintf("Removed %d variants with gnomAD AF > %.1f%%. Remaining: %d variants.\n",
                    sum(exclude_freq), af_threshold * 100, nrow(variants_for_vrs)))
    } else {
        variants_for_vrs <- vcf.df
    }

    if (nrow(variants_for_vrs) == 0) {
        warning("No variants remained after filtering. Skipping VRS query and output generation.")
        return(invisible(list(vus_df = data.frame(), note_df = data.frame())))
    }

    # --- 5. VRS QUERY ---
    
    # Helper: VRS API fetch function
    fetch_vrs_fields <- function(hgvs, retries = 3, timeout_sec = 10) {
      if (is.na(hgvs) || !nzchar(hgvs)) return(c(vrs_id = NA_character_, vrs_seq_id = NA_character_))
      for (attempt in seq_len(retries)) {
        url_str <- httr::modify_url("https://reg.genome.network", path = "vrAllele", query = list(hgvs = hgvs))
        resp <- try(httr::GET(url_str, httr::timeout(timeout_sec)), silent = TRUE)
        if (inherits(resp, "try-error") || httr::http_error(resp)) {
          if (attempt < retries && httr::status_code(resp) %in% c(429, 500:599)) {
            Sys.sleep(0.4 * attempt + stats::runif(1, 0, 0.2))
            next
          } else {
            return(c(vrs_id = NA_character_, vrs_seq_id = NA_character_))
          }
        }
        parsed <- try(jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE), silent = TRUE)
        if (!inherits(parsed, "try-error") && is.list(parsed) && !is.null(parsed[["_id"]])) {
          vrs_id <- parsed[["_id"]]
          vrs_seq_id <- if (!is.null(parsed[["location"]][["sequence_id"]])) parsed[["location"]][["sequence_id"]] else NA_character_
          return(c(vrs_id = vrs_id, vrs_seq_id = vrs_seq_id))
        }
      }
      c(vrs_id = NA_character_, vrs_seq_id = NA_character_)
    }
    
    has_vrs_values <- function(res_vec) {
      if (is.null(res_vec)) return(FALSE)
      any(!is.na(res_vec) & nzchar(unname(res_vec)))
    }
    
    hgvsg_col <- intersect(c("HGVSG", "hgvsg", "koios_hgvsg"), names(variants_for_vrs))[1]
    hgvs_vals <- variants_for_vrs[[hgvsg_col]]
    n <- length(hgvs_vals)
    out <- matrix(NA_character_, nrow = n, ncol = 2, dimnames = list(NULL, c("vrs_id","vrs_seq_id")))

    for (i in seq_len(n)) {
      primary <- hgvs_vals[i]
      res <- fetch_vrs_fields(primary)

      if (!has_vrs_values(res) && "koios_hgvsg" %in% names(variants_for_vrs)) {
        alt <- variants_for_vrs[["koios_hgvsg"]][i]
        if (!is.na(alt) && nzchar(alt) && !identical(alt, primary)) {
          res_alt <- fetch_vrs_fields(alt)
          if (has_vrs_values(res_alt)) res <- res_alt
        }
      }
      out[i, ] <- res
      Sys.sleep(0.05)
    }

    variants_for_vrs$vrs_id <- out[, "vrs_id"]
    variants_for_vrs$vrs_seq_id <- out[, "vrs_seq_id"]

    # --- 6. PHENOPACKETS JSON GENERATION ---

    generate_phenopacket <- function(row_data, clinical_data, vcf_filename) {
        vrs_id <- row_data[["vrs_id"]]
        omop_concept_name <- row_data[["concept_name"]]
        omop_concept_id <- row_data[["concept_id"]]
        hgvsg <- row_data[["hgvsg"]]
        hgvsp_3letter <- row_data[["hgvsp_3letter"]]
        hgvsp_1letter <- row_data[["hgvsp_1letter"]]
        allele_frequency <- row_data[["allele_frequency"]]
        sample_id <- clinical_data$sample_id

        meta_data <- list(
            "phenopacketSchemaVersion" = "2.0",
            "resource" = list(list("id" = "omop_genomic_vocab", "name" = "OMOP Genomic Vocabulary", "namespacePrefix" = "OMOP", "url" = "https://example.com/omop", "version" = "2024", "capability" = "knowledge")),
            "submittedBy" = "KOIOS-VRS-Pipeline"
        )

        subject <- list(
            "id" = sample_id,
            "sex" = list("label" = clinical_data$sex_label, "id" = clinical_data$sex_pato_id),
            "ageOfOnsetToPhenoPacket" = clinical_data$age_of_onset_iso
        )

        # Build extensions list
        extensions <- list(
            list("id" = "omop_classification", "value" = omop_concept_name, "type" = "string", "ontologyClass" = list("id" = omop_concept_id, "label" = omop_concept_name)),
            list("id" = "hgvs_g", "value" = hgvsg, "type" = "string")
        )

        # Add allele frequency if available
        if (!is.null(allele_frequency) && length(allele_frequency) > 0 && !is.na(allele_frequency)) {
            extensions <- c(extensions, list(list("id" = "allele_frequency", "value" = allele_frequency, "type" = "number")))
        }

        # Add protein HGVS annotations if available
        if (!is.na(hgvsp_3letter) && nzchar(hgvsp_3letter)) {
            extensions <- c(extensions, list(list("id" = "hgvs_p_3letter", "value" = hgvsp_3letter, "type" = "string")))
        }

        if (!is.na(hgvsp_1letter) && nzchar(hgvsp_1letter)) {
            extensions <- c(extensions, list(list("id" = "hgvs_p_1letter", "value" = hgvsp_1letter, "type" = "string")))
        }

        variant_interpretation <- list(
            "id" = paste0(sample_id, "_", row_data[["CHROM"]], "_", row_data[["POS"]]),
            "vrsAllele" = list("id" = vrs_id, "type" = "Allele"),
            "extensions" = extensions
        )
        
        interpretations <- list(
            list("id" = "diagnosis_1", "diagnosis" = list("disease" = list("label" = clinical_data$diagnosis_label, "id" = clinical_data$diagnosis_ncit_id))),
            list("id" = "variant_interpretation_1", "variantInterpretation" = variant_interpretation, "progressStatus" = "COMPLETED")
        )
        
        genomic_files <- list(
            list("uri" = vcf_filename, "fileAttributes" = list("fileName" = basename(vcf_filename), "fileFormat" = "VCF", "individualId" = sample_id))
        )
        
        phenopacket <- list(
            "id" = uuid::UUIDgenerate(), 
            "subject" = subject,
            "metaData" = meta_data,
            "interpretations" = interpretations,
            "genomicFiles" = genomic_files
        )
        
        return(phenopacket)
    }

    variants_with_vrs <- variants_for_vrs[!is.na(variants_for_vrs$vrs_id) & nzchar(variants_for_vrs$vrs_id), ]
    phenopacket_list <- lapply(1:nrow(variants_with_vrs), function(i) {
        generate_phenopacket(variants_with_vrs[i, ], clinical_data, vcf_path)
    })

    # Write outputs to output directory
    phenopackets_path <- file.path(output_dir, paste0(output_base_name, "_Phenopackets.json"))
    jsonlite::write_json(phenopacket_list, phenopackets_path, pretty = TRUE, auto_unbox = TRUE)
    cat(sprintf("Phenopackets written to: %s\n", phenopackets_path))

    # --- 7. CSV EXPORT (VUS and Note) ---

    vus_vrs_df <- variants_for_vrs[variants_for_vrs$concept_id == default_vus_id, ]
    note_vrs_df <- variants_for_vrs[variants_for_vrs$concept_id != default_vus_id, ]

    vus_csv_path <- file.path(output_dir, paste0(output_base_name, "_VUS.csv"))
    note_csv_path <- file.path(output_dir, paste0(output_base_name, "_Note.csv"))

    utils::write.csv(vus_vrs_df, file = vus_csv_path, row.names = FALSE, quote = FALSE)
    utils::write.csv(note_vrs_df, file = note_csv_path, row.names = FALSE, quote = FALSE)

    cat(sprintf("VUS CSV written to: %s\n", vus_csv_path))
    cat(sprintf("Note CSV written to: %s\n", note_csv_path))

    # --- 8. GENE-LEVEL MUTATION TABLE ---

    # Extract gene name from concept_name (format: "GENE on GRCh38 chr...")
    extract_gene <- function(concept_name) {
        if (is.na(concept_name) || !nzchar(concept_name)) return(NA_character_)
        # Try to extract gene name before " on GRCh38"
        gene_match <- regmatches(concept_name, regexpr("^[^ ]+", concept_name))
        if (length(gene_match) > 0) return(gene_match[1])
        return(NA_character_)
    }

    variants_for_vrs$gene <- sapply(variants_for_vrs$concept_name, extract_gene)

    # Aggregate mutations by gene
    gene_table <- stats::aggregate(
        cbind(n_mutations = rep(1, nrow(variants_for_vrs))) ~ gene + concept_id,
        data = variants_for_vrs[!is.na(variants_for_vrs$gene), ],
        FUN = length
    )

    # Sort by number of mutations (descending)
    gene_table <- gene_table[order(-gene_table$n_mutations), ]

    gene_table_path <- file.path(output_dir, paste0(output_base_name, "_GeneSummary.csv"))
    utils::write.csv(gene_table, file = gene_table_path, row.names = FALSE, quote = FALSE)
    cat(sprintf("Gene summary table written to: %s\n", gene_table_path))

    # --- 9. FINAL ANNOTATED VCF ---

    # Create final VCF with VRS and KOIOS annotations in INFO field
    final_vcf <- preprocessed$snp_indel

    # Add annotations to INFO field
    for (i in seq_len(nrow(variants_for_vrs))) {
        row <- variants_for_vrs[i, ]

        # Find matching variant in VCF
        vcf_idx <- which(
            final_vcf@fix[, "CHROM"] == row$CHROM &
            final_vcf@fix[, "POS"] == row$POS &
            final_vcf@fix[, "REF"] == row$REF &
            final_vcf@fix[, "ALT"] == row$ALT
        )

        if (length(vcf_idx) > 0) {
            idx <- vcf_idx[1]

            # Build INFO annotations
            info_fields <- character(0)

            if (!is.na(row$vrs_id) && nzchar(row$vrs_id)) {
                info_fields <- c(info_fields, paste0("VRS_ID=", row$vrs_id))
            }

            if (!is.na(row$concept_id) && nzchar(row$concept_id)) {
                info_fields <- c(info_fields, paste0("OMOP_CONCEPT=", row$concept_id))
            }

            if (!is.na(row$gene) && nzchar(row$gene)) {
                info_fields <- c(info_fields, paste0("GENE=", row$gene))
            }

            if (!is.na(row$hgvsg) && nzchar(row$hgvsg)) {
                info_fields <- c(info_fields, paste0("HGVSG=", row$hgvsg))
            }

            if (!is.na(row$hgvsp_3letter) && nzchar(row$hgvsp_3letter)) {
                info_fields <- c(info_fields, paste0("HGVSP_3=", row$hgvsp_3letter))
            }

            if (!is.na(row$hgvsp_1letter) && nzchar(row$hgvsp_1letter)) {
                info_fields <- c(info_fields, paste0("HGVSP_1=", row$hgvsp_1letter))
            }

            # Append to existing INFO or create new
            existing_info <- final_vcf@fix[idx, "INFO"]
            if (length(info_fields) > 0) {
                new_info <- paste(info_fields, collapse = ";")
                if (!is.na(existing_info) && nzchar(existing_info)) {
                    final_vcf@fix[idx, "INFO"] <- paste(existing_info, new_info, sep = ";")
                } else {
                    final_vcf@fix[idx, "INFO"] <- new_info
                }
            }
        }
    }

    # Write final annotated VCF (uncompressed)
    final_vcf_path <- file.path(output_dir, paste0(output_base_name, "_annotated.vcf"))
    temp_final_gz_path <- paste0(final_vcf_path, ".gz")
    vcfR::write.vcf(final_vcf, file = temp_final_gz_path)
    system(sprintf("gunzip -f '%s'", temp_final_gz_path))
    cat(sprintf("Final annotated VCF written to: %s\n", final_vcf_path))

    return(invisible(list(vus_df = vus_vrs_df, note_df = note_vrs_df, gene_table = gene_table)))
}
