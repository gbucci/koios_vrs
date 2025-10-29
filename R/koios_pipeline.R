# koios_pipeline.R

#' KOIOS and VRS Pipeline for Cancer Variants
#'
#' Processes a VCF file, enriches it with KOIOS (HGVSg, ClinGen, OMOP), 
#' queries the VRS API for VRS IDs, filters variants based on frequency, 
#' and generates Phenopackets v2 JSON output along with VUS/Note CSV files.
#'
#' @param vcf_path Path to the input VCF file (must be hg38).
#' @param output_base_name Base name for output files (e.g., 'patient_A' produces 'patient_A_VUS.csv', etc.).
#' @param af_threshold Allele Frequency threshold for filtering (default: 0.01 or 1%). Variants with population AF > threshold are excluded.
#' @param default_vus_id OMOP Concept ID for Variant of Unknown Significance (default: 1028197L).
#' @return A list containing the VUS and Note dataframes, invisibly. Writes output files to disk.
#' @export
koios_vrs_pipeline <- function(vcf_path, output_base_name, af_threshold = 0.01, default_vus_id = 1028197L) {
    
    # --- 1. SETUP AND INPUT VALIDATION ---
    
    clinical_csv_path <- sub("\\.vcf$", ".csv", vcf_path, ignore.case = TRUE)
    if (!file.exists(clinical_csv_path)) {
        stop(sprintf("Clinical metadata file not found: %s.", clinical_csv_path))
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
    
    # --- 2. KOIOS PROCESSING ---
    
    concepts <- KOIOS::loadConcepts()
    vcf <- KOIOS::loadVCF(userVCF = vcf_path)
    ref.df <- KOIOS::loadReference("hg38")
    
    vcf.df <- KOIOS::processVCF(vcf)
    vcf.df <- KOIOS::generateHGVSG(vcf = vcf.df, ref = ref.df)
    vcf.df <- KOIOS::processClinGen(vcf.df, ref = "hg38", progressBar = FALSE)
    vcf.df <- KOIOS::addConcepts(vcf.df, concepts, returnAll = FALSE)
    
    default_vus_name <- "Variant of unknown significance"
    vcf.df$concept_id[is.na(vcf.df$concept_id)] <- default_vus_id
    vcf.df$concept_name[is.na(vcf.df$concept_name) | !nzchar(vcf.df$concept_name)] <- default_vus_name
    
    # --- 3. FILTERING ---
    
    freq_col <- "koios_gnomAD_AF"
    if (freq_col %in% names(vcf.df)) {
        vcf.df[[freq_col]] <- suppressWarnings(as.numeric(vcf.df[[freq_col]]))
        exclude_freq <- !is.na(vcf.df[[freq_col]]) & vcf.df[[freq_col]] > af_threshold
        variants_for_vrs <- vcf.df[!exclude_freq, ]
    } else {
        warning("Column 'koios_gnomAD_AF' not found. Skipping AF filtering.")
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

    jsonlite::write_json(phenopacket_list, paste0(output_base_name, "_Phenopackets.json"), pretty = TRUE, auto_unbox = TRUE)

    # --- 6. CSV EXPORT (VUS and Note) ---

    vus_vrs_df <- variants_for_vrs[variants_for_vrs$concept_id == default_vus_id, ]
    note_vrs_df <- variants_for_vrs[variants_for_vrs$concept_id != default_vus_id, ]

    utils::write.csv(vus_vrs_df, file = paste0(output_base_name, "_VUS.csv"), row.names = FALSE, quote = FALSE)
    utils::write.csv(note_vrs_df, file = paste0(output_base_name, "_Note.csv"), row.names = FALSE, quote = FALSE)
    
    return(invisible(list(vus_df = vus_vrs_df, note_df = note_vrs_df)))
}
