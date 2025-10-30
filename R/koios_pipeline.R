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
#' @return List with 'snp_indel' vcfR object and 'cnv' vcfR object
preprocess_vcf <- function(vcf_obj, min_af = 0.01, min_dp = 100) {

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

    # Apply filters
    pass_af <- is.na(af_vals) | af_vals >= min_af
    pass_dp <- is.na(dp_vals) | dp_vals >= min_dp
    filter_pass <- pass_af & pass_dp

    cat(sprintf("  Removed %d variants with AF < %.1f%%\n", sum(!pass_af), min_af * 100))
    cat(sprintf("  Removed %d variants with DP < %d\n", sum(!pass_dp), min_dp))

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
        # For SNP/INDEL: select allele with highest FAO or AF
        snp_fix <- fix[snp_indel_idx, , drop = FALSE]
        snp_gt <- gt[snp_indel_idx, , drop = FALSE]
        snp_info <- snp_fix[, "INFO"]

        # Extract FAO (Filtered Allele Observation)
        fao_vals <- sapply(snp_info, function(x) {
            fao_match <- regmatches(x, regexec("FAO=([0-9,]+)", x))
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
        cat(sprintf("  Removed %d alternative alleles (kept most representative)\n",
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
#' @param input_dir Input directory containing VCF and clinical CSV files (default: "input").
#' @param output_dir Output directory for all generated files (default: "output").
#' @param default_vus_id OMOP Concept ID for Variant of Unknown Significance (default: 1028197L).
#' @return A list containing the VUS and Note dataframes, invisibly. Writes output files to disk.
#' @export
koios_vrs_pipeline <- function(vcf_path, output_base_name, ref_genome = "hg38", af_threshold = 0.01, min_dp = 100,
                                input_dir = "input", output_dir = "output", default_vus_id = 1028197L) {
    
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

    # Find clinical CSV (same base name as VCF)
    vcf_basename <- basename(vcf_path)
    csv_basename <- sub("\\.vcf$", ".csv", vcf_basename, ignore.case = TRUE)
    clinical_csv_path <- file.path(dirname(vcf_path), csv_basename)

    if (!file.exists(clinical_csv_path)) {
        stop(sprintf("Clinical metadata file not found: %s", clinical_csv_path))
    }
    
    # 1.1. Get Sample ID from VCF
    vcf_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
    sample_names <- colnames(vcf_obj@gt)[-1] 
    if (length(sample_names) != 1) {
        stop("VCF must contain exactly one sample column. Found: ", length(sample_names))
    }
    sample_id <- sample_names[1]
    
    # 1.2. Load Clinical Metadata
    clinical_data <- utils::read.csv(clinical_csv_path, stringsAsFactors = FALSE)
    clinical_data <- clinical_data[clinical_data$sample_id == sample_id, ]
    
    if (nrow(clinical_data) != 1) {
        stop("Clinical metadata not found or does not correspond to the VCF sample ID in the CSV.")
    }

    # --- 1.3. LIFTOVER FROM HG19 TO HG38 IF NEEDED ---

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

        vcfR::write.vcf(vcf_hg38_obj, file = vcf_hg38_path)
        cat(sprintf("Lifted VCF written to: %s\n", vcf_hg38_path))

        # Use the lifted VCF for the rest of the pipeline
        vcf_path <- vcf_hg38_path
        # Reload vcf_obj with lifted coordinates
        vcf_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
    }

    # --- 1.4. VCF PREPROCESSING ---

    # Preprocess VCF: filter quality, remove multi-allelic, separate CNV
    preprocessed <- preprocess_vcf(vcf_obj, min_af = af_threshold, min_dp = min_dp)

    if (is.null(preprocessed$snp_indel) || nrow(preprocessed$snp_indel@fix) == 0) {
        warning("No SNP/INDEL variants remained after preprocessing. Skipping KOIOS processing.")
        return(invisible(list(vus_df = data.frame(), note_df = data.frame())))
    }

    # Write preprocessed VCF for KOIOS in output directory
    preprocessed_vcf_filename <- sub("\\.vcf$", "_preprocessed.vcf", basename(vcf_path), ignore.case = TRUE)
    preprocessed_vcf_path <- file.path(output_dir, preprocessed_vcf_filename)
    vcfR::write.vcf(preprocessed$snp_indel, file = preprocessed_vcf_path)
    cat(sprintf("Preprocessed VCF written to: %s\n\n", preprocessed_vcf_path))

    # Optional: Write CNV to separate file for future processing
    if (!is.null(preprocessed$cnv) && nrow(preprocessed$cnv@fix) > 0) {
        cnv_vcf_filename <- sub("\\.vcf$", "_CNV.vcf", basename(vcf_path), ignore.case = TRUE)
        cnv_vcf_path <- file.path(output_dir, cnv_vcf_filename)
        vcfR::write.vcf(preprocessed$cnv, file = cnv_vcf_path)
        cat(sprintf("CNV variants written to: %s\n\n", cnv_vcf_path))
    }

    # --- 2. KOIOS PROCESSING ---

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

    # --- 3. POST-KOIOS FILTERING ---
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

    # --- 4. VRS QUERY ---
    
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

    # --- 5. PHENOPACKETS JSON GENERATION ---

    generate_phenopacket <- function(row_data, clinical_data, vcf_filename) {
        vrs_id <- row_data[["vrs_id"]]
        omop_concept_name <- row_data[["concept_name"]]
        omop_concept_id <- row_data[["concept_id"]]
        hgvsg <- row_data[["hgvsg"]]
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
        
        variant_interpretation <- list(
            "id" = paste0(sample_id, "_", row_data[["CHROM"]], "_", row_data[["POS"]]),
            "vrsAllele" = list("id" = vrs_id, "type" = "Allele"),
            "extensions" = list(
                list("id" = "omop_classification", "value" = omop_concept_name, "type" = "string", "ontologyClass" = list("id" = omop_concept_id, "label" = omop_concept_name)),
                list("id" = "hgvs_g", "value" = hgvsg, "type" = "string")
            )
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

    # --- 6. CSV EXPORT (VUS and Note) ---

    vus_vrs_df <- variants_for_vrs[variants_for_vrs$concept_id == default_vus_id, ]
    note_vrs_df <- variants_for_vrs[variants_for_vrs$concept_id != default_vus_id, ]

    vus_csv_path <- file.path(output_dir, paste0(output_base_name, "_VUS.csv"))
    note_csv_path <- file.path(output_dir, paste0(output_base_name, "_Note.csv"))

    utils::write.csv(vus_vrs_df, file = vus_csv_path, row.names = FALSE, quote = FALSE)
    utils::write.csv(note_vrs_df, file = note_csv_path, row.names = FALSE, quote = FALSE)

    cat(sprintf("VUS CSV written to: %s\n", vus_csv_path))
    cat(sprintf("Note CSV written to: %s\n", note_csv_path))
    
    return(invisible(list(vus_df = vus_vrs_df, note_df = note_vrs_df)))
}
