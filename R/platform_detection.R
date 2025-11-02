# platform_detection.R
# Functions to detect sequencing platform and adjust processing parameters

#' Detect Sequencing Platform from VCF Header
#'
#' Detects the sequencing platform (Ion Torrent, Illumina, etc.) from VCF metadata
#'
#' @param vcf_path Path to VCF file
#' @return List with platform name, caller, version, and recommended parameters
#' @export
detect_sequencing_platform <- function(vcf_path) {

    # Read VCF header
    header_lines <- readLines(vcf_path, n = 200)
    header_lines <- header_lines[grepl("^##", header_lines)]

    platform_info <- list(
        platform = "unknown",
        caller = "unknown",
        version = NA_character_,
        panel = NA_character_,
        params = list()
    )

    # --- ION TORRENT DETECTION ---
    # Pattern 1: ##source="tvc 5.12-27 (e1619b46ee) - Torrent Variant Caller"
    tvc_match <- grep("source.*tvc.*Torrent Variant Caller", header_lines, value = TRUE, ignore.case = TRUE)
    if (length(tvc_match) > 0) {
        platform_info$platform <- "IonTorrent"
        platform_info$caller <- "TVC"

        # Extract version
        version_match <- regmatches(tvc_match[1], regexpr("tvc [0-9.-]+", tvc_match[1], ignore.case = TRUE))
        if (length(version_match) > 0) {
            platform_info$version <- sub("tvc ", "", version_match[1], ignore.case = TRUE)
        }

        # Detect panel from header
        panel_match <- grep("TargetRegions|TargetPanel|##Panel", header_lines, value = TRUE, ignore.case = TRUE)
        if (length(panel_match) > 0) {
            # Try to extract Oncomine panel name
            if (grepl("Oncomine", panel_match[1], ignore.case = TRUE)) {
                onco_match <- regmatches(panel_match[1], regexpr("Oncomine[^,;>\"]+", panel_match[1], ignore.case = TRUE))
                if (length(onco_match) > 0) {
                    platform_info$panel <- onco_match[1]
                }
            }
        }

        # Ion Torrent specific parameters
        platform_info$params <- list(
            min_af = 0.01,      # Ion Torrent has higher background noise
            min_dp = 100,       # Good coverage typically needed
            min_fao = 10,       # Filter out low-quality variants
            use_fao = TRUE,     # FAO field is specific to Ion Torrent
            multiallelic_common = TRUE  # Ion Torrent often reports multi-allelic sites
        )

        return(platform_info)
    }

    # --- ILLUMINA DETECTION ---
    # Pattern 1: ##source=strelka, ##source=GATK, ##source=VarDict, etc.
    illumina_callers <- c("strelka", "gatk", "mutect", "varscan", "vardict", "bcftools", "freebayes", "lofreq")
    for (caller in illumina_callers) {
        caller_match <- grep(paste0("source.*", caller), header_lines, value = TRUE, ignore.case = TRUE)
        if (length(caller_match) > 0) {
            platform_info$platform <- "Illumina"
            platform_info$caller <- toupper(caller)

            # Extract version if available
            version_match <- regmatches(caller_match[1], regexpr("[0-9]+\\.[0-9]+[.0-9]*", caller_match[1]))
            if (length(version_match) > 0) {
                platform_info$version <- version_match[1]
            }

            # Check for Myriapod panel
            if (grepl("Myriapod|Myriad", paste(header_lines, collapse = " "), ignore.case = TRUE)) {
                platform_info$panel <- "Myriapod"
            }

            # Illumina specific parameters
            platform_info$params <- list(
                min_af = 0.05,      # Illumina typically has lower noise
                min_dp = 50,        # Lower DP threshold acceptable
                min_fao = NA,       # FAO not typically used
                use_fao = FALSE,
                use_ad = TRUE,      # Use AD (Allelic Depth) instead
                multiallelic_common = FALSE
            )

            return(platform_info)
        }
    }

    # Pattern 2: Check instrument/platform fields
    platform_match <- grep("##platform=|##instrument=", header_lines, value = TRUE, ignore.case = TRUE)
    if (length(platform_match) > 0) {
        if (grepl("illumina|nextseq|novaseq|miseq|hiseq", platform_match[1], ignore.case = TRUE)) {
            platform_info$platform <- "Illumina"
            platform_info$params <- list(
                min_af = 0.05,
                min_dp = 50,
                min_fao = NA,
                use_fao = FALSE,
                use_ad = TRUE,
                multiallelic_common = FALSE
            )
            return(platform_info)
        }

        if (grepl("ion|torrent|proton|s5|genestudio", platform_match[1], ignore.case = TRUE)) {
            platform_info$platform <- "IonTorrent"
            platform_info$params <- list(
                min_af = 0.01,
                min_dp = 100,
                min_fao = 10,
                use_fao = TRUE,
                multiallelic_common = TRUE
            )
            return(platform_info)
        }
    }

    # --- UNKNOWN PLATFORM ---
    # Use conservative defaults
    warning("Could not detect sequencing platform. Using conservative default parameters.")
    platform_info$params <- list(
        min_af = 0.05,      # More stringent
        min_dp = 100,       # Higher coverage
        min_fao = 10,       # Conservative
        use_fao = FALSE,
        multiallelic_common = FALSE
    )

    return(platform_info)
}

#' Get Platform-Specific Preprocessing Parameters
#'
#' Returns optimized parameters for VCF preprocessing based on platform
#'
#' @param platform_info Output from detect_sequencing_platform()
#' @param override_params Optional list to override specific parameters
#' @return List of preprocessing parameters
#' @export
get_platform_parameters <- function(platform_info, override_params = NULL) {

    params <- platform_info$params

    # Override with user-specified parameters
    if (!is.null(override_params)) {
        for (param_name in names(override_params)) {
            params[[param_name]] <- override_params[[param_name]]
        }
    }

    return(params)
}

#' Print Platform Detection Summary
#'
#' Prints a summary of detected platform and parameters
#'
#' @param platform_info Output from detect_sequencing_platform()
#' @export
print_platform_summary <- function(platform_info) {
    cat("\n=== SEQUENCING PLATFORM DETECTION ===\n")
    cat(sprintf("Platform: %s\n", platform_info$platform))
    cat(sprintf("Variant Caller: %s\n", platform_info$caller))

    if (!is.na(platform_info$version)) {
        cat(sprintf("Version: %s\n", platform_info$version))
    }

    if (!is.na(platform_info$panel)) {
        cat(sprintf("Panel: %s\n", platform_info$panel))
    }

    cat("\nRecommended Parameters:\n")
    cat(sprintf("  min_af: %.2f%%\n", platform_info$params$min_af * 100))
    cat(sprintf("  min_dp: %d\n", platform_info$params$min_dp))

    if (!is.na(platform_info$params$min_fao)) {
        cat(sprintf("  min_fao: %d\n", platform_info$params$min_fao))
    }

    cat("\n")
}
