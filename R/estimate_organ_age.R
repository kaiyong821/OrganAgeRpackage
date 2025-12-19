#' Estimate Organ Ages from Proteomic Data
#'
#' Predicts biological age for multiple organs using SomaLogic proteomic data
#'
#' @param data Data frame with protein expression data containing:
#'   \itemize{
#'     \item Age column (chronological age)
#'     \item Sex column (0 = male, 1 = female)
#'     \item ID column (sample identifier)
#'     \item Protein columns starting with "seq." (RFU units, NOT log-transformed)
#'   }
#' @param age_col Character. Name of age column. Default: "AGE"
#' @param sex_col Character. Name of sex column (0=male, 1=female). Default: "Sex_F"
#' @param id_col Character. Name of ID column. Default: "HELPFul_ID"
#' @param protein_prefix Character. Protein column prefix. Default: "seq."
#' @param lowess_frac Numeric. LOWESS smoothing span (0-1). Default: 2/3
#' @param lowess_iter Integer. LOWESS iterations. Default: 5
#' @param verbose Logical. Print progress messages. Default: TRUE
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item ID: Sample identifier
#'     \item Age: Chronological age
#'     \item Sex_F: Sex (0=male, 1=female)
#'     \item Predicted_[Organ]_Age: Predicted organ age
#'     \item AgeGap_[Organ]: Age gap (predicted - cohort LOWESS expected)
#'     \item AgeGap_zscore_[Organ]: Z-scored age gap (cohort-normalized)
#'   }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Log10-transforms protein expression data
#'   \item Z-score normalizes using KADRC healthy control parameters
#'   \item Predicts organ age using 500 bootstrap LASSO models
#'   \item Fits cohort-specific LOWESS curve (predicted age ~ chronological age)
#'   \item Calculates age gaps: predicted age - LOWESS expected age
#'   \item Z-scores age gaps within each organ using cohort distribution
#' }
#'
#' Note: Cognition-related organs are excluded from output.
#'
#' @importFrom dplyr %>% filter mutate select rename across all_of
#' @importFrom stats lowess approx
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' \dontrun{
#' result <- estimate_organ_age(data = my_proteomics_data)
#' }
estimate_organ_age <- function(data,
                               age_col = "AGE",
                               sex_col = "Sex_F",
                               id_col = "HELPFul_ID",
                               protein_prefix = "seq.",
                               lowess_frac = 2/3,
                               lowess_iter = 5,
                               verbose = TRUE) {

  # Suppress R CMD check notes
  ID <- Age <- Sex_F <- organ <- protein_seqid <- seq_col <- NULL

  # Important data format notice
  if (verbose) {
    cat("\n===============================================\n")
    cat("OrganAgeR: Organ Age Prediction\n")
    cat("===============================================\n\n")
    cat("IMPORTANT: Please ensure your SomaScan data:\n")
    cat("  - Is from SomaScan 7K assay\n")
    cat("  - Is exported from anmlSMP.ADat file\n")
    cat("  - Contains RAW RFU values (NOT log-transformed)\n")
    cat("  - If using other versions, please lift to 7K first\n\n")
  }

  # Input validation
  required_cols <- c(age_col, sex_col, id_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  unique_sex <- unique(data[[sex_col]])
  if (!all(unique_sex %in% c(0, 1, NA))) {
    stop(sex_col, " must contain only 0 (male), 1 (female), or NA")
  }

  seq_cols <- grep(paste0("^", protein_prefix), names(data), value = TRUE)
  if (length(seq_cols) == 0) {
    stop("No protein columns found with prefix: ", protein_prefix)
  }

  # Load model parameters
  scaler_path <- system.file("extdata", "KADRC_protein_scaler_params.csv",
                             package = "OrganAgeR")

  if (!file.exists(scaler_path)) {
    stop("Scaler file not found. Please ensure OrganAgeR is properly installed.")
  }

  kadrc_scalers <- utils::read.csv(scaler_path, stringsAsFactors = FALSE)

  # Load model coefficients from multiple split files
  extdata_path <- system.file("extdata", package = "OrganAgeR")
  model_files <- list.files(extdata_path,
                            pattern = "^models_.*\\.csv$",
                            full.names = TRUE)

  if (length(model_files) == 0) {
    stop("Model files not found. Please ensure OrganAgeR is properly installed.")
  }

  if (verbose) {
    cat("Loading", length(model_files), "organ model files...\n")
  }

  model_coefficients <- do.call(rbind, lapply(model_files, function(f) {
    utils::read.csv(f, stringsAsFactors = FALSE)
  }))

  # Exclude Cognition organs
  cognition_organs <- grep("^Cognition", unique(model_coefficients$organ), value = TRUE)
  if (length(cognition_organs) > 0) {
    model_coefficients <- model_coefficients[!model_coefficients$organ %in% cognition_organs, ]
  }

  kadrc_scalers <- kadrc_scalers[!kadrc_scalers$organ %in% cognition_organs, ]

  # Data preprocessing
  data_processed <- data %>%
    dplyr::rename(
      ID = !!rlang::sym(id_col),
      Age = !!rlang::sym(age_col),
      Sex_F = !!rlang::sym(sex_col)
    ) %>%
    dplyr::select(ID, Age, Sex_F, dplyr::all_of(seq_cols))

  data_processed <- data_processed %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(seq_cols), log10))

  # Initialize results
  result_df <- data.frame(
    ID = data_processed$ID,
    Age = data_processed$Age,
    Sex_F = data_processed$Sex_F,
    stringsAsFactors = FALSE
  )

  # Predict organ ages
  organs <- unique(model_coefficients$organ)

  if (verbose) {
    cat("Starting prediction for", length(organs), "organs...\n\n")
  }

  for (organ_name in organs) {

    if (verbose) cat("Computing", organ_name, "age...")

    organ_scalers <- kadrc_scalers %>%
      dplyr::filter(organ == organ_name) %>%
      dplyr::mutate(seq_col = paste0(protein_prefix, gsub("-", ".", protein_seqid)))

    if (nrow(organ_scalers) == 0) {
      if (verbose) cat(" SKIPPED (no scaler parameters)\n")
      next
    }

    available_proteins <- intersect(organ_scalers$seq_col, seq_cols)

    # Get actual number of proteins used in models (exclude all-zero coefficients)
    organ_models_temp <- model_coefficients %>%
      dplyr::filter(organ == organ_name)

    # Identify features with non-zero coefficients
    python_cols_temp <- names(organ_models_temp)[!names(organ_models_temp) %in%
                                                   c("organ", "bootstrap_seed", "y_intercept")]
    coef_temp <- as.matrix(organ_models_temp[, python_cols_temp])

    # Use loop to check each column for non-zero coefficients
    non_zero_mask <- logical(length(python_cols_temp))
    for (i in seq_along(python_cols_temp)) {
      col_vals <- as.numeric(coef_temp[, i])
      non_zero_mask[i] <- any(abs(col_vals) > 1e-10, na.rm = TRUE)
    }

    non_zero_features <- python_cols_temp[non_zero_mask]

    # Convert to R feature names
    non_zero_r_names <- non_zero_features
    for (i in seq_along(non_zero_r_names)) {
      if (non_zero_r_names[i] != "Sex_F") {
        clean_name <- sub("^X", "", non_zero_r_names[i])
        non_zero_r_names[i] <- paste0(protein_prefix, clean_name)
      }
    }

    # Count active proteins (excluding Sex_F)
    active_protein_features <- non_zero_r_names[non_zero_r_names != "Sex_F"]
    required_proteins <- length(active_protein_features)

    # Check how many active proteins are actually available in the data
    available_active_proteins <- intersect(active_protein_features, seq_cols)
    n_available <- length(available_active_proteins)

    # Safety check
    if (is.na(required_proteins) || is.null(required_proteins) || required_proteins == 0) {
      if (verbose) cat(" SKIPPED (invalid protein count)\n")
      next
    }

    # Check if all active proteins are available
    if (n_available < required_proteins) {
      if (verbose) {
        cat(" SKIPPED (only ", n_available, "/",
            required_proteins, " proteins available)\n", sep = "")
      }
      next
    } else {
      # Show protein availability for all organs
      if (verbose) {
        cat(" (", n_available, "/", required_proteins,
            " proteins available)", sep = "")
      }
    }

    X_organ <- as.matrix(data_processed[, available_proteins, drop = FALSE])

    kadrc_lookup <- organ_scalers %>%
      dplyr::filter(seq_col %in% available_proteins)
    kadrc_lookup <- kadrc_lookup[match(available_proteins, kadrc_lookup$seq_col), ]

    X_organ_scaled <- sweep(X_organ, 2, kadrc_lookup$kadrc_mean, "-")
    X_organ_scaled <- sweep(X_organ_scaled, 2, kadrc_lookup$kadrc_std, "/")

    X_with_sex <- cbind(Sex_F = data_processed$Sex_F, X_organ_scaled)

    organ_models <- model_coefficients %>%
      dplyr::filter(organ == organ_name)

    intercepts <- as.numeric(gsub("\\[|\\]", "", organ_models$y_intercept))

    python_feature_cols <- names(organ_models)
    python_feature_cols <- python_feature_cols[!python_feature_cols %in%
                                                 c("organ", "bootstrap_seed", "y_intercept")]

    r_feature_names <- python_feature_cols
    for (i in seq_along(r_feature_names)) {
      if (r_feature_names[i] != "Sex_F") {
        clean_name <- sub("^X", "", r_feature_names[i])
        r_feature_names[i] <- paste0(protein_prefix, clean_name)
      }
    }

    coef_matrix <- as.matrix(organ_models[, python_feature_cols])
    colnames(coef_matrix) <- r_feature_names

    # Exclude features with all-zero coefficients from prediction (with floating point tolerance)
    active_features_mask <- apply(coef_matrix, 2, function(x) any(abs(x) > 1e-10))
    active_features <- colnames(coef_matrix)[active_features_mask]

    model_features <- r_feature_names[r_feature_names %in% colnames(X_with_sex) &
                                        r_feature_names %in% active_features]

    if (length(model_features) == 0) {
      if (verbose) cat(" SKIPPED (no matched features)\n")
      next
    }

    X_for_pred <- X_with_sex[, model_features, drop = FALSE]
    coef_matrix_subset <- coef_matrix[, model_features, drop = FALSE]

    bootstrap_predictions <- X_for_pred %*% t(coef_matrix_subset)
    bootstrap_predictions <- sweep(bootstrap_predictions, 2, intercepts, "+")

    mean_predictions <- rowMeans(bootstrap_predictions)

    pred_col <- paste0("Predicted_", organ_name, "_Age")
    result_df[[pred_col]] <- mean_predictions

    lowess_fit <- stats::lowess(
      x = result_df$Age,
      y = mean_predictions,
      f = lowess_frac,
      iter = lowess_iter
    )

    expected_age <- suppressWarnings(
      stats::approx(
        x = lowess_fit$x,
        y = lowess_fit$y,
        xout = result_df$Age,
        rule = 2,
        ties = mean
      )$y
    )

    age_gap <- mean_predictions - expected_age
    age_gap_mean <- mean(age_gap, na.rm = TRUE)
    age_gap_sd <- sqrt(mean((age_gap - age_gap_mean)^2, na.rm = TRUE))
    age_gap_zscore <- (age_gap - age_gap_mean) / age_gap_sd

    result_df[[paste0("AgeGap_", organ_name)]] <- age_gap
    result_df[[paste0("AgeGap_zscore_", organ_name)]] <- age_gap_zscore

    if (verbose) cat(" DONE\n")
  }

  n_organs <- sum(grepl("^Predicted_", names(result_df)))

  if (verbose) {
    cat("\n===============================================\n")
    cat("Successfully predicted", n_organs, "organs\n")
    cat("===============================================\n\n")
  }

  return(result_df)
}
